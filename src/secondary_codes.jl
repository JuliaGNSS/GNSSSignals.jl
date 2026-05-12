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

Receivers should query the secondary code via [`get_secondary_code`](@ref),
which returns the appropriate `SecondaryCode` instance for the signal.

The internal API for evaluating the secondary at a given secondary-chip
index is `secondary_value`. It is typed-dispatched so the
no-secondary case folds to the multiplicative identity at compile time.
"""
abstract type SecondaryCode end

"""
    NoSecondaryCode <: SecondaryCode

Marker type for signals without a secondary code. All
`secondary_value` lookups return `1`, which the inner code-generation
loops multiply by — the compiler elides the multiplication entirely.
"""
struct NoSecondaryCode <: SecondaryCode end

@inline secondary_code_length(::NoSecondaryCode) = 1
# Returning `true` (not `1`) so that `x * true === x` preserves `eltype` for
# every numeric element type — no widening to `Int64` on a code matrix of
# `Int16` or `Float32`.
@inline secondary_value(::NoSecondaryCode, ::Integer, ::Integer) = true

"""
    SharedSecondaryCode{N, T<:Integer} <: SecondaryCode

A short secondary code that is identical across all PRNs. Stored as an
`NTuple{N, T}` so the values are inlined into the type and indexed by the
compiler at zero runtime cost.

`code` holds the secondary chip values as a tuple of `±1` (or another small
integer alphabet, depending on the signal's spec).
"""
struct SharedSecondaryCode{N, T<:Integer} <: SecondaryCode
    code::NTuple{N, T}

    # Inner constructor that takes the chips as varargs (at least one).
    # Defining any inner constructor suppresses Julia's auto-generated
    # outer constructor; the `c1::T, crest::T...` signature pins `T`
    # from the first chip and counts `N` from the splat, keeping both
    # type parameters bound for Aqua's `unbound_args` check.
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

A long secondary code whose values differ per PRN. Stored as a
`(length, num_prns)` matrix — same layout convention as the primary code
matrix on the signal.

Used by GPS L1C (1800-chip overlay) and signals with similar per-PRN
overlays.
"""
struct PerPRNSecondaryCode{T<:Integer, M<:AbstractMatrix{T}} <: SecondaryCode
    codes::M
end

@inline secondary_code_length(s::PerPRNSecondaryCode) = size(s.codes, 1)
@inline function secondary_value(s::PerPRNSecondaryCode, prn::Integer, secondary_index::Integer)
    @inbounds s.codes[mod(secondary_index, size(s.codes, 1)) + 1, prn]
end
