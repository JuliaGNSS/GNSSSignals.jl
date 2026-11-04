# `GNSSSignals.CodeLUT` — internal submodule, not part of the public API; use `gen_code!` /
# `gen_code` / `code_engine`. Its names (`BPSK`, `BOC`, `TMBOC`, `CBOC`, …) live inside it so they
# do not clash with GNSSSignals' own modulations. (Plain comment, not a docstring, so Documenter's
# checkdocs doesn't require this internal module in the manual.)
module CodeLUT

using SIMD

# ---- backends ----
abstract type Backend end
struct AVX512   <: Backend end   # vpermb over a 64-chip sliding window (W = 64)
struct AVX2     <: Backend end   # vpshufb over two independent 16-chip windows (W = 32)
struct Neon     <: Backend end   # tbl1 over a single 16-chip window (W = 16, AArch64)
struct Portable <: Backend end   # scalar fallback (any CPU)

backend_name(::AVX512)   = "AVX-512"
backend_name(::AVX2)     = "AVX2"
backend_name(::Neon)     = "NEON"
backend_name(::Portable) = "portable"

# Window width / SIMD lane count per backend.
_vwidth(::AVX512)   = Val(64)
_vwidth(::AVX2)     = Val(32)
_vwidth(::Neon)     = Val(16)
_vwidth(::Portable) = Val(1)

"""
    CodeTable(chips::AbstractVector{<:Integer})

Holds a GNSS spreading code of length `L` as `Int8` chips (values are stored verbatim;
pass ±1 for a standard correlation replica). Internally keeps a copy padded with its own
first 63 chips so any 64-chip window `chips[base : base+63]` is a single contiguous load.
"""
struct CodeTable{V<:AbstractVector{Int8}}
    chips::V    # length L
    padded::V   # length L + WINDOW_PAD
    length::Int
end

# vpermb reads a 64-chip window; the last valid base is L-1, reading up to index L+62
# (0-based). Pad by 63 so that load is always in-bounds.
const WINDOW_PAD = 63

# Zero-copy view tables: the auto-generated outer constructor `CodeTable(chips, padded, length)`
# infers `V`, so unit-stride `SubArray{Int8}` column views (e.g. of `SignalLUT.padded`) give a
# view-backed table with no copy; the SIMD `VecRange{W}` loads lower to the same contiguous vload
# as for a `Vector{Int8}`. Passing `Vector{Int8}` (the bake path) gives an owning table.

function CodeTable(chips::AbstractVector{<:Integer})
    L = length(chips)
    L > 0 || throw(ArgumentError("code length must be positive"))
    c = Int8.(chips)
    padded = vcat(c, c[1:min(WINDOW_PAD, L)])
    # if L < WINDOW_PAD, repeat until we have L + WINDOW_PAD entries
    while length(padded) < L + WINDOW_PAD
        padded = vcat(padded, c[1:min(L + WINDOW_PAD - length(padded), L)])
    end
    CodeTable(c, padded, L)
end

Base.length(t::CodeTable) = t.length

include("permute.jl")
include("kernel.jl")
include("iterate.jl")
include("modulation.jl")
include("generator.jl")

# ---- backend selection ----
@static if Sys.ARCH in (:x86_64, :i686)
    # The AVX-512 branch is gated on the PRECOMPILE-time `HOST_FEATURES` const: the
    # `vpermb`/`Vec{64,Int8}` `FillEngine512` path only compiles (LLVM aborts legalising `<64 x i8>`
    # ops) when the pkgimage was precompiled on an AVX-512 host, so `default_backend()` must never
    # return `AVX512()` unless the const enabled it.
    #
    # The const may still OVER-claim: the pkgimage can be loaded on a CPU without VBMI
    # (JULIA_CPU_TARGET=generic, shared/NFS depot, Docker/CI) where `vpermb` SIGILLs (#104). So
    # inside that branch we consult the RUNTIME feature set (`RUNTIME_FEATURES`, filled in
    # `__init__`) and demote AVX512 -> AVX2/Portable. The AVX2-only branch likewise re-checks avx2
    # at runtime and demotes to Portable (#126).
    @static if HOST_FEATURES.avx512vbmi
        # Delegates to `_select_backend` (the unit-tested policy) on the RUNTIME features so the
        # shipping path can't drift from the tested helper (#129). One Ref load per call.
        default_backend() = _select_backend(RUNTIME_FEATURES[])
    elseif HOST_FEATURES.avx2
        # The AVX-512 path is not compiled here, so never return AVX512(). AVX2 can over-claim
        # like #104: with JULIA_CPU_TARGET=generic the AVX2 `_pshufb` (an llvmcall forcing +avx2)
        # SIGILLs on a non-AVX2 CPU (#126), so re-check at runtime.
        default_backend() = RUNTIME_FEATURES[].avx2 ? AVX2() : Portable()
    else
        default_backend() = Portable()
    end

    # Re-validate CPU features on the running machine: `HOST_FEATURES` was baked at precompile and
    # may over-claim on a relocated pkgimage (see permute.jl / #104).
    __init__() = _refresh_host_features!()
elseif Sys.ARCH === :aarch64
    default_backend() = Neon()
else
    default_backend() = Portable()
end
# Length-aware: AVX2/NEON widen the phase vector to Int32 for tables > typemax(Int16) (slower
# than Int16, but ~20× over scalar); only longer tables (unreachable for GNSS) fall to Portable.
function default_backend(table::CodeTable)
    be = default_backend()
    (be isa Union{AVX2,Neon} && table.length > typemax(Int32)) ? Portable() : be
end

end # module CodeLUT
