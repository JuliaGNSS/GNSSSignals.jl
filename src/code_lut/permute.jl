# Low-level permute intrinsic + CPU-feature detection.
#
# The whole speed-up rests on one instruction: AVX-512 VBMI `vpermb` treats a 512-bit
# register as a 64-entry Int8 lookup table and gathers all 64 output lanes in a single
# op (`out[i] = table[index[i] & 63]`). We use it to read a 64-chip *window* of the PRN
# code at relative indices 0..63 (see kernel.jl for why a 64-entry window suffices).
#
# Only the Int8 `vpermb` is wired up here — that is all the code resampler needs. The
# AVX2 (`vpshufb`) and NEON (`tbl4`) windowed paths are future work; see default_backend.

@static if Sys.ARCH in (:x86_64, :i686)

# ---- CPU feature detection (cpuid leaf 7 + xgetbv), mirrors SinCosLUT ----
@inline function _cpuid(leaf::UInt32, subleaf::UInt32)
    Base.llvmcall(
        """
        %s = call { i32, i32, i32, i32 } asm sideeffect "cpuid",
             "={ax},={bx},={cx},={dx},{ax},{cx},~{dirflag},~{fpsr},~{flags}"(i32 %0, i32 %1)
        %a = extractvalue { i32, i32, i32, i32 } %s, 0
        %b = extractvalue { i32, i32, i32, i32 } %s, 1
        %c = extractvalue { i32, i32, i32, i32 } %s, 2
        %d = extractvalue { i32, i32, i32, i32 } %s, 3
        %r0 = insertvalue [4 x i32] undef, i32 %a, 0
        %r1 = insertvalue [4 x i32] %r0, i32 %b, 1
        %r2 = insertvalue [4 x i32] %r1, i32 %c, 2
        %r3 = insertvalue [4 x i32] %r2, i32 %d, 3
        ret [4 x i32] %r3
        """,
        Tuple{UInt32,UInt32,UInt32,UInt32}, Tuple{UInt32,UInt32}, leaf, subleaf)
end
@inline _bit(x::UInt32, n) = (x >> n) & 0x1 == 0x1
@inline _xcr0() = Base.llvmcall(
    """
    %r = call i32 asm sideeffect "xgetbv", "={ax},{cx},~{dx},~{dirflag},~{fpsr},~{flags}"(i32 %0)
    ret i32 %r
    """, UInt32, Tuple{UInt32}, UInt32(0))

function _x86_features()
    _, _, ecx1, _ = _cpuid(UInt32(1), UInt32(0))
    osxsave = _bit(ecx1, 27)
    xcr0 = osxsave ? _xcr0() : UInt32(0)
    avx_os    = osxsave && _bit(xcr0, 1) && _bit(xcr0, 2)
    avx512_os = avx_os  && _bit(xcr0, 5) && _bit(xcr0, 6) && _bit(xcr0, 7)
    _, ebx, ecx, _ = _cpuid(UInt32(7), UInt32(0))
    (avx2       = avx_os    && _bit(ebx, 5),
     avx512vbmi = avx512_os && _bit(ecx, 1))
end

# Detected once at precompile and baked into a const (like VectorizationBase.jl). Kept only as
# a seed/hint — it can OVER-claim: if the pkgimage is precompiled on an AVX-512-VBMI host and
# later loaded on a CPU without VBMI (JULIA_CPU_TARGET=generic, shared/NFS depot, Docker/CI),
# pkgimage validation checks only cloned code targets, not serialized const data, so this stale
# const would still say `avx512vbmi = true`. Backend selection therefore consults the RUNTIME
# feature set below, not this const, so we never emit `vpermb` on a non-VBMI CPU (SIGILL, #104).
const HOST_FEATURES = _x86_features()

# The CPU features of the machine ACTUALLY running, re-detected once per session in `__init__`
# (see `default_backend`). Seeded with the precompile-time const so a read before `__init__`
# (unusual) is still valid; `__init__` overwrites it with a live `cpuid` on the running CPU.
const RUNTIME_FEATURES = Ref(HOST_FEATURES)

# Re-run CPU-feature detection on the executing machine and store it. Called from `__init__`.
_refresh_host_features!() = (RUNTIME_FEATURES[] = _x86_features(); nothing)

# Pure backend selection from a CPU-feature set (a plain NamedTuple), split out so it is
# unit-testable without a live CPU: a non-VBMI set must NEVER yield `AVX512()` (whose `vpermb`
# would SIGILL), while a VBMI-capable set may (keeping the fast path for capable CPUs). #104.
_select_backend(features) =
    features.avx512vbmi ? AVX512() :
    features.avx2       ? AVX2()   : Portable()

# ---- vpermb: 64-entry Int8 in-register table lookup ----
# `alwaysinline` so llvmcall emits the bare permute instead of a real call + spill.
const _PERMB_IR = """
declare <64 x i8> @llvm.x86.avx512.permvar.qi.512(<64 x i8>, <64 x i8>)
define <64 x i8> @entry(<64 x i8> %0, <64 x i8> %1) #0 {
  %r = call <64 x i8> @llvm.x86.avx512.permvar.qi.512(<64 x i8> %0, <64 x i8> %1)
  ret <64 x i8> %r }
attributes #0 = { alwaysinline "target-features"="+avx512vbmi,+avx512bw,+avx512f" }
"""

# `out[i] = table[index[i] & 63]`. Indices outside 0..63 wrap (low 6 bits); our windowed
# relative indices are always in 0..63, so no masking is needed.
@inline _permute(table::Vec{64,Int8}, index::Vec{64,Int8}) =
    Vec{64,Int8}(Base.llvmcall((_PERMB_IR, "entry"), NTuple{64,VecElement{Int8}},
        Tuple{NTuple{64,VecElement{Int8}},NTuple{64,VecElement{Int8}}}, table.data, index.data))

# ---- vpshufb: per-128-bit-lane 16-entry Int8 table lookup (AVX2) ----
# Shuffles each 128-bit lane independently using the low 4 bits of each index (bit 7
# zeroes the output lane). We exploit the lane independence to look up two separate
# 16-chip windows in one instruction; our indices are 0..15, so no zeroing occurs.
const _PSHUFB_IR = """
declare <32 x i8> @llvm.x86.avx2.pshuf.b(<32 x i8>, <32 x i8>)
define <32 x i8> @entry(<32 x i8> %0, <32 x i8> %1) #0 {
  %r = call <32 x i8> @llvm.x86.avx2.pshuf.b(<32 x i8> %0, <32 x i8> %1)
  ret <32 x i8> %r }
attributes #0 = { alwaysinline "target-features"="+avx2" }
"""
@inline _pshufb(table::Vec{32,Int8}, index::Vec{32,Int8}) =
    Vec{32,Int8}(Base.llvmcall((_PSHUFB_IR, "entry"), NTuple{32,VecElement{Int8}},
        Tuple{NTuple{32,VecElement{Int8}},NTuple{32,VecElement{Int8}}}, table.data, index.data))

end # @static x86

# ---- tbl1: 16-entry Int8 in-register table lookup (AArch64 NEON) ----
# `tbl1` treats a 128-bit register as a 16-entry Int8 lookup table and gathers all 16
# output lanes (`out[i] = table[index[i]]`; out-of-range indices ≥ 16 return 0). It is the
# single-window analogue of one of AVX2's two `vpshufb` halves: 16 chips → 16 lanes.
#
# NOT @static-guarded: defining this llvmcall on x86 is harmless — it only compiles when
# CALLED, and the Neon backend is never selected/called on x86 (see default_backend). It is
# validated on macOS-ARM CI. Mirrors the CI-proven SinCosLUT `tbl4` shape (permute_neon.jl).
const _TBL1_IR = """
declare <16 x i8> @llvm.aarch64.neon.tbl1(<16 x i8>, <16 x i8>)
define <16 x i8> @entry(<16 x i8> %t, <16 x i8> %i) #0 {
  %r = call <16 x i8> @llvm.aarch64.neon.tbl1(<16 x i8> %t, <16 x i8> %i)
  ret <16 x i8> %r }
attributes #0 = { alwaysinline "target-features"="+neon" }
"""
@inline _tbl1(tbl::Vec{16,Int8}, index::Vec{16,Int8}) =
    Vec{16,Int8}(Base.llvmcall((_TBL1_IR, "entry"), NTuple{16,VecElement{Int8}},
        Tuple{NTuple{16,VecElement{Int8}},NTuple{16,VecElement{Int8}}}, tbl.data, index.data))
