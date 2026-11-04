# The chip a phase falls in: the phase counts chips for every modulation.
get_floored_phase(::Modulation, phase) = floor(Int, phase)

"""
$(SIGNATURES)

Get the code value at a given phase for a specific PRN.

Returns the spreading code value (including subcarrier modulation for BOC signals)
at the specified code phase. The phase is automatically wrapped to the code length.

# Arguments
- `phase`: Code phase in chips

Returns a value of type [`get_code_type`](@ref)`(signal)`: `Int`, or `Float32` for CBOC.
The stored `Int8` chip is widened, so correlations and sums over `get_code` values cannot
overflow.

# Examples
```julia-repl
julia> get_code(GPSL1CA(), 0.0, 1)
1
julia> get_code(GPSL1CA(), 1200.3, 1)
-1
julia> get_code.(GPSL1CA(), 0:1022, 1)  # Full code period
```
"""
function get_code(signal::T, phase, prn::Integer) where {T<:AbstractGNSSSignal}
    get_code(get_modulation(T), signal, phase, prn)
end

# An integer chip is widened to `Int`: the built-in matrices store `Int8`, whose sums and
# products wrap around after 127.
@inline _widen_chip(chip::Integer) = Int(chip)
@inline _widen_chip(chip) = chip

# Multiply the primary chip by its secondary-code chip. The absolute chip is
# either wrapped by get_code or supplied in range by get_code_unsafe.
Base.@propagate_inbounds function _primary_secondary_code(signal, absolute_chip, prn)
    primary_length = size(signal.codes, 1)
    sec = get_secondary_code(signal)
    chip_idx = mod(absolute_chip, primary_length)
    sec_idx = div(absolute_chip, primary_length)
    _widen_chip(get_code_at_index(signal, chip_idx, prn)) * secondary_value(sec, prn, sec_idx)
end

@inline _modulate_code(::BPSK, code, phase) = code
@inline _modulate_code(modulation::BOC, code, phase) =
    code * get_subcarrier_code(modulation, phase)

function get_code(modulation::Union{BOC,BPSK}, signal::AbstractGNSSSignal, phase, prn::Integer)
    floored_phase = get_floored_phase(modulation, phase)
    primary_length = size(signal.codes, 1)
    sec_len = secondary_code_length(get_secondary_code(signal))
    absolute_chip = mod(floored_phase, primary_length * sec_len)
    _modulate_code(modulation, _primary_secondary_code(signal, absolute_chip, prn), phase)
end

Base.@propagate_inbounds function get_code_at_index(
    signal::AbstractGNSSSignal,
    phase::Integer,
    prn::Integer,
)
    signal.codes[phase+1, prn]
end

"""
$(SIGNATURES)

Get the code value at a given phase without bounds checking.

This is a faster version of [`get_code`](@ref) that skips the modulo operation.
The phase must be within `[0, code_length * secondary_code_length)`.

!!! warning
    Using phases outside the valid range results in undefined behavior.

# Examples
```julia-repl
julia> get_code_unsafe(GPSL1CA(), 500.0, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(
    signal::S,
    phase,
    prn::Integer,
) where {S<:AbstractGNSSSignal}
    get_code_unsafe(get_modulation(S), signal, phase, prn)
end

# The caller supplies a phase within the primary × secondary code length.
Base.@propagate_inbounds function get_code_unsafe(
    modulation::Union{BOC,BPSK},
    signal::AbstractGNSSSignal,
    phase,
    prn::Integer,
)
    floored_phase = get_floored_phase(modulation, phase)
    _modulate_code(modulation, _primary_secondary_code(signal, floored_phase, prn), phase)
end

"""
$(SIGNATURES)

Get the full code matrix for a GNSS signal.

Returns the primary codes as a matrix where each column represents a PRN.
Built-in signals store `Int8` chips; sampled generation reads a separate baked LUT.

!!! warning "Widen before correlating"
    `Int8` arithmetic wraps around silently, and a product or sum of two `Int8`s stays
    `Int8`. `dot(a, b)`, `a' * b`, `codes' * codes` and `reduce(+, a)` over more than 127
    chips therefore return wrong values without an error (`sum` is safe: it widens).
    Widen first, e.g. `Int.(get_codes(signal))`, or use [`get_code`](@ref), which
    returns `Int`.

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

Read codes of the given type from a binary file, returning a `Matrix{type}` of size
`(code_length, num_prns)`.

# Examples
```julia-repl
julia> read_in_codes(Int8, "/data/gpsl1cacodes.bin", 32, 1023)
```
"""
function read_in_codes(::Type{T}, filename, num_prns, code_length) where {T}
    # Two-argument `open` and explicit `close`, not the `do`-block form: its keyword
    # plumbing is the one call Julia's `--trim=safe` verifier cannot resolve.
    file_stream = open(filename, "r")
    try
        read!(file_stream, Array{T}(undef, code_length, num_prns))
    finally
        close(file_stream)
    end
end
