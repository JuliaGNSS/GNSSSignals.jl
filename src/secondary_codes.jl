"""
    SecondaryCode

Abstract supertype for a signal's secondary (overlay) code. A secondary code
modulates the primary code at one chip per primary code period.

Three concrete forms cover all currently-defined GNSS signals:

- [`NoSecondaryCode`](@ref) — signals without a secondary code (GPS L1 C/A,
  Galileo E1B, etc.). Evaluation is a compile-time no-op.
- [`SharedSecondaryCode`](@ref) — short overlay that is the same for all
  PRNs (GPS L5 NH10/NH20, Galileo E1 CS25).
- [`PerPRNSecondaryCode`](@ref) — long per-PRN overlay (GPS L1C overlay,
  1800 chips per PRN from a per-PRN LFSR).

Query it via [`get_secondary_code`](@ref) and read its chips with
[`secondary_value`](@ref GNSSSignals.secondary_value). It is type-dispatched so the
no-secondary case folds to the multiplicative identity at compile time.
"""
abstract type SecondaryCode end

"""
    secondary_value(secondary_code::SecondaryCode, prn, index) -> ±1

The chip of `secondary_code` at the 0-based secondary-chip `index` for satellite
`prn`: the sign the primary code period `index` (counted from the start of the
secondary code) is overlaid with. `index` wraps modulo the code's length, and
`prn` is read only by a [`PerPRNSecondaryCode`](@ref). The value has the code's
element type, and is `true` for a [`NoSecondaryCode`](@ref), so that multiplying
by it keeps any numeric element type.

```julia-repl
julia> sec = get_secondary_code(GPSL5I());

julia> [GNSSSignals.secondary_value(sec, 1, k) for k in 0:9]'
1×10 adjoint(::Vector{Int8}) with eltype Int8:
 1  1  1  1  -1  -1  1  -1  1  -1
```
"""
function secondary_value end

"""
    NoSecondaryCode <: SecondaryCode

Marker type for signals without a secondary code. All `secondary_value` lookups return
`1`, so the compiler elides the multiplication.
"""
struct NoSecondaryCode <: SecondaryCode end

@inline secondary_code_length(::NoSecondaryCode) = 1
# `true` (not `1`) so that `x * true === x` preserves `eltype` for every numeric type.
@inline secondary_value(::NoSecondaryCode, ::Integer, ::Integer) = true

"""
    SharedSecondaryCode{N, T<:Integer} <: SecondaryCode

A short secondary code that is identical across all PRNs. `code` holds the chips as an
`NTuple{N, T}` of `±1` (or another small integer alphabet, per the signal's spec), so they
are inlined into the type.
"""
struct SharedSecondaryCode{N, T<:Integer} <: SecondaryCode
    code::NTuple{N, T}

    # Varargs constructor (at least one chip): `c1::T, crest::T...` pins `T` and counts `N`,
    # keeping both type parameters bound for Aqua's `unbound_args` check.
    function SharedSecondaryCode(c1::T, crest::T...) where {T<:Integer}
        N = 1 + length(crest)
        return new{N, T}((c1, crest...))
    end
end

@inline secondary_code_length(::SharedSecondaryCode{N}) where {N} = N
@inline function secondary_value(s::SharedSecondaryCode{N}, ::Integer, secondary_index::Integer) where {N}
    @inbounds s.code[mod(secondary_index, N) + 1]
end

"""
    PerPRNSecondaryCode{T<:Integer} <: SecondaryCode

A long secondary code whose values differ per PRN, stored as a `(length, num_prns)`
matrix like the primary code matrix. Used by GPS L1C (1800-chip overlay).
"""
struct PerPRNSecondaryCode{T<:Integer, M<:AbstractMatrix{T}} <: SecondaryCode
    codes::M
end

@inline secondary_code_length(s::PerPRNSecondaryCode) = size(s.codes, 1)
@inline function secondary_value(s::PerPRNSecondaryCode, prn::Integer, secondary_index::Integer)
    @inbounds s.codes[mod(secondary_index, size(s.codes, 1)) + 1, prn]
end
