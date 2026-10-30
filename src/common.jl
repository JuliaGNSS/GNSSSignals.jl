"""
$(SIGNATURES)

Get the full code matrix for a GNSS signal.

Returns the codes as a matrix where each column represents a PRN.

# Arguments
- `signal`: A GNSS signal instance (e.g., `GPSL1CA()`, `GPSL5I()`, `GalileoE1B()`)

# Returns
- `Matrix`: Code matrix of size `(code_length, num_prns)`

# Examples
```julia-repl
julia> codes = get_codes(GPSL1CA())
julia> size(codes)
(1023, 37)
```
"""
function get_codes(signal::AbstractGNSSSignal)
    signal.codes
end

"""
$(SIGNATURES)

Widen a primary-code matrix from its on-disk / LFSR-generated `Int8`
representation to `Int16`.

Chip values are ±1 and would fit in `Int8`, but on x86_64 / AVX2 hardware
storing as `Int16` is materially faster for `gen_code!`: the inner store
loop emits a clean `vpbroadcastw` + `vmovq` pattern, while Int8 storage
triggers an `shl 8 + or` byte-packing antipattern (3 extra μops per chip)
when the buffer is also `Int8`. Storing Int8 chips into an Int16 buffer
recovers the clean codegen but still loses ~14 % because the `movsx`
load chain runs slower than `movzx`.

Probed alternatives that did **not** recover Int8 perf on AVX2:
- `@simd ivdep` annotation on the inner store loop (no measurable
  effect — LLVM has already made its codegen choice for `Val`-known
  fixed-trip loops)
- `unsafe_store!` of a replicated `UInt32` to bypass LLVM's
  pattern-matcher (helped Int8/Int8 by ~5 % but still ~12 % slower
  than Int16/Int16, and a regression for Int16/Int16)
- `SIMD.jl` `Vec{N,T}` broadcast-store (slower at every NUM_INNER
  measured because the abstraction forces an xmm materialization for
  small N)

Memory cost of widening is small in absolute terms — the largest
current matrix (GPS L5-I, 10230 × 37) is 757 KB at Int16. The trade-off
may invert on AVX-512 or non-x86 hardware; revisit if you have access
to such platforms.
"""
widen_codes_to_storage(codes::AbstractMatrix) = Int16.(codes)

"""
$(SIGNATURES)

Get the `Symbol` identifier of a GNSS signal — the machine-readable per-signal
key, e.g. `:GPSL1CA`, `:GalileoE1B`, `:GalileoE5aI`.

This is the finest identity level: the same PRN measured on two bands of one
constellation (e.g. GPS L1 and L5) has two distinct signal ids.

Defaults to `nameof` of the signal type; override `get_signal_id(::Type{MySignal})`
to pin a different symbol. The BOC(1,1) approximation types keep their own ids
(`:GalileoE1B_BOC11`), distinct from the signals they approximate. Works on an instance or a type and folds to a
compile-time constant.

Distinct from [`get_signal_name`](@ref), which is a human-readable display
string (`"GPS L1 C/A"`), and coarser than [`get_band_id`](@ref) (`:L1`).

```julia-repl
julia> get_signal_id(GPSL1CA())
:GPSL1CA

julia> get_signal_id(GalileoE1B)
:GalileoE1B
```
"""
@inline get_signal_id(::Type{S}) where {S<:AbstractGNSSSignal} = nameof(S)
@inline get_signal_id(s::AbstractGNSSSignal) = get_signal_id(typeof(s))

"""
$(SIGNATURES)

Get the human-readable name of a GNSS signal, e.g. `"GPS L1 C/A"`,
`"Galileo E5a-Q"`.

This is the display counterpart to [`get_signal_id`](@ref), at the same (finest)
identity level: it names one specific transmission, component included. The
spelling follows the ICD rather than the type name, so it is meant for log lines,
plot labels and user-facing output — prefer the `Symbol` id when branching or
keying. The BOC(1,1) approximation types name themselves as such (`"Galileo E1B
(BOC(1,1) approximation)"`), so a label never silently claims to be the exact
signal.

Genuinely per signal, so each concrete signal defines one method on its own type
(`get_signal_name(::Type{<:GPSL1CA}) = "GPS L1 C/A"`, in that signal's file),
mirroring [`get_band`](@ref); the instance method forwards to the type, so the
name is available without constructing a signal.

Finer than [`get_band_name`](@ref) (`"L1"`) and [`get_constellation_name`](@ref)
(`"GPS"`).

```julia-repl
julia> get_signal_name(GPSL1CA())
"GPS L1 C/A"

julia> get_signal_name(GalileoE5aQ)
"Galileo E5a-Q"
```
"""
function get_signal_name end

@inline get_signal_name(s::AbstractGNSSSignal) = get_signal_name(typeof(s))

"""
$(SIGNATURES)

Get the `Symbol` identifier of the GNSS constellation a signal belongs to —
`:GPS`, `:Galileo` or `:BeiDou`.

This is the coarsest identity level: every signal of one constellation, across
all bands and PRNs, maps to the same id (GPS L1 C/A, GPS L5I and GPS L1C are all
`:GPS`). It's the natural key for logging, dictionaries and constellation-level
branching; for a type-level test prefer `signal isa AbstractGPSSignal`.

Defined once per constellation on the abstract signal type
(`get_constellation_id(::Type{<:AbstractGPSSignal}) = :GPS`), so every concrete
signal inherits it without constructing a value. Works on an instance or a type
and folds to a compile-time constant.

Distinct from [`get_time_system`](@ref) (the constellation's reference time
scale, `GPST()`/`GST()`/`BDT()`), and coarser than [`get_band_id`](@ref) (`:L1`)
and [`get_signal_id`](@ref) (`:GPSL1CA`).

```julia-repl
julia> get_constellation_id(GPSL1CA())
:GPS

julia> get_constellation_id(GalileoE1B)
:Galileo
```
"""
@inline get_constellation_id(::Type{<:AbstractGPSSignal}) = :GPS
@inline get_constellation_id(::Type{<:AbstractGalileoSignal}) = :Galileo
@inline get_constellation_id(::Type{<:AbstractBeiDouSignal}) = :BeiDou
@inline get_constellation_id(s::AbstractGNSSSignal) = get_constellation_id(typeof(s))

"""
$(SIGNATURES)

Get the human-readable name of the GNSS constellation a signal belongs to —
`"GPS"`, `"Galileo"` or `"BeiDou"`.

The display counterpart to [`get_constellation_id`](@ref) — same granularity
(every signal of one constellation returns the same string), but a `String`
meant for log lines and user-facing output rather than a key to branch or
dictionary on (prefer the `Symbol` id for that; comparing symbols is cheaper
and the name is free to be respelled).

Stated as a literal once per constellation, on the abstract signal type
(`get_constellation_name(::Type{<:AbstractGPSSignal}) = "GPS"`) as
[`get_constellation_id`](@ref) is, so every concrete signal inherits the name
without constructing a value, the call allocates nothing and it folds to a
compile-time constant. Each literal still spells exactly what `String` of the id
yields, and a test pins that, so the two layers cannot drift apart.

A constellation added later needs no method of its own: the
`String(get_constellation_id(S))` fallback names it for free — at the cost of
allocating a fresh `String` on every call, for the reason spelled out at
[`get_band_name`](@ref), whose derived fallback is built the same way. State a
literal, as the constellations here do, if you are naming one in a hot loop, or
wherever the display name differs from the id anyway: a future `:NavIC` reading
`"NavIC (IRNSS)"`.

Works on an instance or a type.

Coarser than [`get_band_name`](@ref) (`"L1"`) and [`get_signal_name`](@ref)
(`"GPS L1 C/A"`).

```julia-repl
julia> get_constellation_name(GPSL1CA())
"GPS"

julia> get_constellation_name(GalileoE1B)
"Galileo"
```
"""
@inline get_constellation_name(::Type{<:AbstractGPSSignal}) = "GPS"
@inline get_constellation_name(::Type{<:AbstractGalileoSignal}) = "Galileo"
@inline get_constellation_name(::Type{<:AbstractBeiDouSignal}) = "BeiDou"
@inline get_constellation_name(::Type{S}) where {S<:AbstractGNSSSignal} =
    String(get_constellation_id(S))
@inline get_constellation_name(s::AbstractGNSSSignal) = get_constellation_name(typeof(s))

# Each concrete signal defines these on `::Type{<:Signal}` (per-signal files);
# these forward an instance to its type.
@inline get_code_length(s::AbstractGNSSSignal) = get_code_length(typeof(s))

"""
    get_chip_rate(signal) -> Frequency

The nominal chip rate of a signal's primary ranging code: the chips per second the
satellite broadcasts, as its ICD states it (1.023 MHz for GPS L1 C/A, 10.23 MHz for GPS
L5). It is a property of the signal; the rate a receiver's replica runs at is its *code
frequency*, this rate shifted by the code Doppler, which [`gen_code!`](@ref) and
[`GNSSSignals.code_engine`](@ref) take as their `code_frequency` argument (defaulting to
this one). Defined per concrete signal type; works on an instance or a type.
"""
@inline get_chip_rate(s::AbstractGNSSSignal) = get_chip_rate(typeof(s))

"""
    get_bit_rate(signal) -> Frequency
    get_bit_rate(signal, prn) -> Frequency

The bit rate of a signal's navigation message as tracked: the rate of the navigation
bits a tracking loop demodulates, one sign each. On a signal with forward error
correction that is the symbol rate of the code, not the bit rate of the message it
carries: GPS L5-I is tracked at 100 bit/s while its CNAV message carries 50 bit/s. Where
the message is not FEC-coded (GPS L1 C/A's LNAV, BeiDou D1) the two coincide. `0 Hz` for
a dataless pilot.

With a `prn`, the rate that satellite broadcasts. It differs from the signal's only on
BeiDou B1I and B3I, whose GEO satellites ([`is_beidou_geo`](@ref GNSSSignals.is_beidou_geo))
broadcast the D2 message at 500 Hz where the signal's rate is D1's 50 Hz. Pass the PRN
wherever the satellite is known.

Defined per concrete signal type; works on an instance or a type.

```julia-repl
julia> get_bit_rate(BeiDouB1I(), 6), get_bit_rate(BeiDouB1I(), 1)
(50 Hz, 500 Hz)
```
"""
@inline get_bit_rate(s::AbstractGNSSSignal) = get_bit_rate(typeof(s))
@inline get_bit_rate(::Type{S}, ::Integer) where {S<:AbstractGNSSSignal} =
    get_bit_rate(S)
@inline get_bit_rate(s::AbstractGNSSSignal, prn::Integer) =
    get_bit_rate(typeof(s), prn)

"""
    get_code_blocks_per_bit(signal) -> Int
    get_code_blocks_per_bit(signal, prn) -> Int

The number of primary code periods (code blocks) a navigation bit spans, as tracked
([`get_bit_rate`](@ref)): `get_chip_rate / (get_code_length · get_bit_rate)`, a whole
number for every data signal — 20 for GPS L1 C/A and Galileo E5a-I, 10 for GPS L5-I, 5 for
BeiDou B2a-I, 4 for Galileo E5b-I, 1 where a bit spans one code period (Galileo E1B, GPS
L1C-D, …). With a `prn`, that satellite's ([`get_bit_rate`](@ref)`(signal, prn)`): 2 on a
BeiDou B1I or B3I GEO satellite, whose D2 bits are ten times shorter.

`0` for a dataless pilot, which broadcasts no bits, as its bit rate is `0 Hz`: `iszero`
tells a pilot apart on either, and a caller that needs no bit constraint for a pilot gets
none (`max(get_secondary_code_length(signal), blocks)`). Works on an instance or a type.

```julia-repl
julia> get_code_blocks_per_bit(GPSL1CA())
20

julia> get_code_blocks_per_bit(BeiDouB1I(), 6), get_code_blocks_per_bit(BeiDouB1I(), 1)
(20, 2)
```
"""
@inline get_code_blocks_per_bit(::Type{S}) where {S<:AbstractGNSSSignal} =
    _code_blocks_per_bit(S, get_bit_rate(S))
@inline get_code_blocks_per_bit(::Type{S}, prn::Integer) where {S<:AbstractGNSSSignal} =
    _code_blocks_per_bit(S, get_bit_rate(S, prn))
@inline get_code_blocks_per_bit(s::AbstractGNSSSignal) =
    get_code_blocks_per_bit(typeof(s))
@inline get_code_blocks_per_bit(s::AbstractGNSSSignal, prn::Integer) =
    get_code_blocks_per_bit(typeof(s), prn)

# In integer Hz, so that a bit rate that does not divide the chip rate into whole blocks
# fails loudly rather than rounding. A rate may be stated as a float, such as n · 1.023 MHz
# for a non-integer n; it is converted to an `Int` first, and one that is no whole number of
# Hz is an error.
function _code_blocks_per_bit(::Type{S}, bit_rate) where {S<:AbstractGNSSSignal}
    iszero(bit_rate) && return 0
    blocks, remainder = divrem(
        _whole_hertz(S, "chip rate", get_chip_rate(S)),
        get_code_length(S) * _whole_hertz(S, "bit rate", bit_rate),
    )
    iszero(remainder) || throw(
        ArgumentError(get_signal_name(S) * "'s bits span no whole number of code blocks"),
    )
    blocks
end

function _whole_hertz(::Type{S}, what, rate) where {S<:AbstractGNSSSignal}
    hertz = ustrip(Hz, rate)
    isinteger(hertz) || throw(ArgumentError(
        get_signal_name(S) * "'s " * what * " of " * string(hertz) * " Hz is no whole number of Hz"))
    Int(hertz)
end

"""
$(SIGNATURES)

Carrier phase, in **cycles**, of the signal's component relative to its band's
in-phase (I) carrier reference — the carrier the ICD designates as the "in-phase"
component of that band (GPS L1/L2: the P(Y) carrier; GPS L5: the I5 carrier;
Galileo E5a: the E5a-I carrier). Two signals on one band share this reference, so
their offsets are directly comparable. A quarter cycle is 90°; `cispi(2 * offset)`
rotates by the offset, and `2π * offset` gives it in radians.

Several bands multiplex components in phase quadrature (QPSK), one on the in-phase
carrier and one on the quadrature carrier. That 90° offset lives purely in the
carrier — the spreading codes are real-valued — so it is not otherwise recoverable
from the signal definition. A receiver that tracks such components jointly off one
shared carrier loop needs this to keep each on its own decision axis (e.g.
de-rotating the data prompt when the loop is locked to the pilot); without it a
quadrature component collapses onto the orthogonal axis and never demodulates.

The sign follows each ICD's own convention, so it differs by constellation:

- **GPS** puts the civil components on the *quadrature* carrier, *lagging* the
  band's P(Y) in-phase reference by 90° → `−0.25`: `GPSL1CA` (IS-GPS-200N §3.3.1.5.1),
  `GPSL2CM`/`GPSL2CL` (IS-GPS-200N §3.3.1.5.1; nominal — CNAV Type 10 bit 273 can
  command L2C in-phase, `0.0`), and `GPSL5Q` (IS-GPS-705 §3.3.1.5). `GPSL5I` is the
  L5 in-phase reference (`0.0`).
- **GPS L1C** (`GPSL1C_D`, `GPSL1C_P`) rides the *same* P(Y) in-phase carrier
  (IS-GPS-800J §3.2.1.6.1) → `0.0`; hence L1C sits 90° off C/A on the same band.
- **Galileo E5aQ** *leads* the E5a-I reference by 90° → `+0.25` (OS SIS ICD Eq. 1,
  `I·cos − Q·sin`); `GalileoE5aI` is the reference (`0.0`).
- **Galileo E5bQ** likewise *leads* its own sideband's `GalileoE5bI` reference by 90°
  → `+0.25`: the E5 AltBOC feeds each sideband as `(I + jQ)` (OS SIS ICD Eq. 2), and
  the ICD blesses processing E5a and E5b as two independent QPSK signals (§2.3.1.2).
  `GalileoE5aQP` shares the E5a carrier but the ICD states no phase relation for it
  (§2.3.1.4), so it reports the `0.0` default — no offset is claimed, not one of `0`.
- **Galileo E1B/E1C** are both on the E1 in-phase carrier → `0.0`; their relative
  180° anti-phase is carried in the CBOC code (`get_modulation`), not the carrier.
- **Galileo E6B/E6C** likewise share one carrier component, but here the 180°
  *is* the carrier's: OS SIS ICD Eq. 10 subtracts the pilot,
  `(e_E6-B − e_E6-C)/√2`, and the plain-BPSK E6-C code carries no sign of its own
  (E1C's does, in `CBOC`'s `boc2_sign`), so `GalileoE6C` is `0.5` against the
  `GalileoE6B` reference (`0.0`).
- **BeiDou B2a** follows the same complex-envelope convention as Galileo
  (BDS-SIS-ICD-B2a-1.0 Eq. 4-2/4-3, `s_data + j·s_pilot` with `I·cos − Q·sin`), so
  the `BeiDouB2aQ` pilot *leads* its `BeiDouB2aI` reference by 90° → `+0.25`.
- **BeiDou B1C** is likewise `s_data + j·s_pilot` (BDS-SIS-ICD-B1C-1.0 Eq. 4-3):
  the ICD's own phase table (Table 4-2) puts the data at 0° and the pilot's
  BOC(1,1) arm at 90°, so `BeiDouB1C_P` — whose replica is that BOC(1,1) arm —
  is `+0.25` against the `BeiDouB1C_D` reference (`0.0`).

Default `0.0` covers every single-component signal and in-phase component above
(including `BeiDouB1I`, `BeiDouB3I` and `BeiDouB2bI`, each its band's I component
per its ICD).
"""
@inline get_carrier_phase_offset(::Type{<:AbstractGNSSSignal}) = 0.0
@inline get_carrier_phase_offset(s::AbstractGNSSSignal) = get_carrier_phase_offset(typeof(s))

"""
$(SIGNATURES)

Nominal power of a signal component, as a **dimensionless linear power**: its share
of its composite signal, with composites on the same constellation and band scaled
against one another. The value is the ICD power split itself, so the components of
one composite sum to `1`.

A loop tracking several components of one satellite weights each discriminator
output by its component's power — `w = map(get_relative_power, signals)`,
normalised. Only ratios matter, so elevation, satellite block, free-space loss and
antenna gain cancel; the absolute ICD figures are worst-case floors several dB
below what a receiver sees, and are deliberately not exposed.

| constellation, band | composite (the unit, `1`) | components |
|:---|:---|:---|
| GPS L1 | L1C | `GPSL1C_P` `0.75`, `GPSL1C_D` `0.25`; `GPSL1CA` `0.708` |
| GPS L2 | L2C | `GPSL2CM` `0.5`, `GPSL2CL` `0.5` |
| GPS L5 | L5 | `GPSL5I` `0.5`, `GPSL5Q` `0.5` |
| Galileo E1 | E1 | `GalileoE1B` `0.5`, `GalileoE1C` `0.5` |
| Galileo E5a | E5a | `GalileoE5aI` `0.5`, `GalileoE5aQ` `0.5`; `GalileoE5aQP` `0.282` |
| Galileo E5b | E5b | `GalileoE5bI` `0.5`, `GalileoE5bQ` `0.5` |
| Galileo E6 | E6 | `GalileoE6B` `0.5`, `GalileoE6C` `0.5` |
| BeiDou B1C | B1C | `BeiDouB1C_D` `0.25`, `BeiDouB1C_P` `0.75` |
| BeiDou B2a | B2a | `BeiDouB2aI` `0.5`, `BeiDouB2aQ` `0.5` |
| BeiDou B1I / B3I / B2b | each its own | `BeiDouB1I` / `BeiDouB3I` / `BeiDouB2bI`, `1.0` each |

Every composite splits evenly except the two modern L1-carrier ones, which are
75/25 pilot/data: GPS L1C (IS-GPS-800J, Table 3.2-1) and BeiDou B1C
(BDS-SIS-ICD-B1C-1.0 §4.2.2, data:pilot 1:3). These are payload design constants,
fixed independently of satellite block and elevation, so each is stated as the
plain split and comes out exact. `BeiDouB1I`, `BeiDouB3I` and `BeiDouB2bI` are
single-component: each is BeiDou's only open-service component on its carrier (the
quadrature counterparts are not in the OS ICDs), so each is its own composite,
`1.0`. That the Galileo E5b components share `BeiDouB2bI`'s carrier does not enter
into it — the scale is per constellation and band, so a BeiDou value is only
comparable with another BeiDou value.

`GPSL1CA` and `GalileoE5aQP` are the two values that are not bare splits, and the
reason the scale spans a band rather than a composite: each is a *separate* signal on
a carrier its composite already occupies, not a further component of that composite,
so each is scaled against it. C/A is `−158.5 dBW` (IS-GPS-200N, Table 3-Va) against
L1C's `−157.0 dBW`, i.e. `10^(-0.15)`; E5a-QP is `−160.75 dBW` against E5a's
`−155.25 dBW` (Galileo OS SIS ICD v2.2, Table 13), i.e. `10^(-0.55)`. That lets a
joint L1 or E5a loop weight all of a carrier's signals against each other directly.

Values from different rows are **not** comparable: a ratio across constellations or
bands would have to account for front-end gain, filter loss and noise figure, which
can swamp a 1–2 dB difference in transmitted power. `GalileoE1B` and `GPSL5I` both
being `0.5` says only that each takes half of its own composite. Weight across bands
from measured C/N₀ instead.

```julia-repl
julia> get_relative_power(GPSL1C_P()), get_relative_power(GPSL1C_D())
(0.75, 0.25)

julia> get_relative_power(GPSL1CA())
0.7079457843841379
```

Defined per concrete signal type beside [`get_band`](@ref), with no generic default,
so a signal without a value raises a `MethodError`. Works on an instance or a type
and folds to a compile-time constant.
"""
@inline get_relative_power(s::AbstractGNSSSignal) = get_relative_power(typeof(s))

"""
$(SIGNATURES)

Generate a sampled code signal for a given PRN.

Allocates a new `Vector{Int8}` and returns the spreading code sampled at the specified
sampling frequency. Thin allocating wrapper over the in-place [`gen_code!`](@ref); see it for
the Int8 output semantics (±1, or the CBOC integer approximation), fractional-phase support,
and the `sampling_frequency ≥ code_frequency · subchip_factor` requirement.

# Arguments
- `num_samples`: Number of samples to generate
- `signal`: GNSS signal instance (e.g., `GPSL1CA()`, `GPSL5I()`, `GalileoE1B()`)
- `prn`: PRN number of the satellite
- `sampling_frequency`: Sampling frequency (must be ≥ code frequency · subchip_factor)
- `code_frequency`: Code frequency the replica runs at (default: the signal's nominal chip
  rate, [`get_chip_rate`](@ref))
- `start_phase`: Initial code phase in chips (default: 0.0)
- `start_index`: Index offset (default: 0)

# Returns
- `Vector{Int8}`: Sampled code signal

# Examples
```julia-repl
julia> using Unitful: MHz
julia> sampled_code = gen_code(4000, GPSL1CA(), 1, 4MHz)
julia> length(sampled_code)
4000
```
"""
function gen_code(
    num_samples::Integer,
    signal::AbstractGNSSSignal,
    prn::Integer,
    sampling_frequency::Frequency,
    code_frequency::Frequency = get_chip_rate(signal),
    start_phase = 0.0,
    start_index::Integer = 0,
)
    code = zeros(Int8, num_samples)
    gen_code!(code, signal, prn, sampling_frequency, code_frequency, start_phase, start_index)
end

"""
$(SIGNATURES)

Get the ratio of chip rate to center frequency.

This ratio is used to compute the code Doppler from the carrier Doppler.

# Arguments
- `signal`: A GNSS signal instance

# Returns
- `Float64`: The code-to-center frequency ratio

# Examples
```julia-repl
julia> get_chip_rate_to_center_frequency_ratio(GPSL1CA())
0.0006493506493506494
```
"""
@inline function get_chip_rate_to_center_frequency_ratio(signal::AbstractGNSSSignal)
    get_chip_rate(signal) / get_center_frequency(signal)
end

"""
$(SIGNATURES)

Get the length of the secondary code.

# Arguments
- `signal`: A GNSS signal instance

# Returns
- `Int`: Secondary code length (1 if no secondary code)

# Examples
```julia-repl
julia> get_secondary_code_length(GPSL1CA())
1
julia> get_secondary_code_length(GPSL5I())
10
```
"""
@inline function get_secondary_code_length(signal::AbstractGNSSSignal)
    secondary_code_length(get_secondary_code(signal))
end
