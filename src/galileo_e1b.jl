"""
    GalileoE1B{C} <: AbstractGNSS{C}

Galileo E1B signal type.

Galileo E1B uses CBOC(6,1,1/11) modulation with a 4092 chip code at 1.023 Mcps,
transmitted on the E1 carrier frequency of 1575.42 MHz.
"""
struct GalileoE1B{C<:AbstractMatrix} <: AbstractGNSS{C}
    codes::C
end

get_modulation(::Type{<:GalileoE1B}) = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11)
@inline get_modulation(::GalileoE1B) = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11)

get_system_string(s::GalileoE1B) = "GalileoE1B"

function read_from_documentation(raw_code)
    raw_code_without_spaces = replace(replace(raw_code, " " => ""), "\n" => "")
    code_hex_array = map(x -> parse(UInt16, x; base = 16), collect(raw_code_without_spaces))
    code_bit_string = string(string.(code_hex_array, base = 2, pad = 4)...)
    map(x -> parse(Int8, x; base = 2), collect(code_bit_string)) .* Int8(2) .- Int8(1)
end

function read_galileo_e1b_codes()
    read_in_codes(
        Int8,
        joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_galileo_e1b.bin"),
        50,
        4092,
    )
end

function GalileoE1B()
    GalileoE1B(Int16.(read_galileo_e1b_codes()))
end

"""
$(SIGNATURES)

Get the code length for Galileo E1B.

# Returns
- `Int`: 4092 chips

# Examples
```julia-repl
julia> get_code_length(GalileoE1B())
4092
```
"""
@inline function get_code_length(galileo_e1b::GalileoE1B)
    4092
end

"""
$(SIGNATURES)

Get the secondary code for Galileo E1B.

Galileo E1B has no secondary code, so this returns 1.

# Returns
- `Int`: 1 (no secondary code)

# Examples
```julia-repl
julia> get_secondary_code(GalileoE1B())
1
```
"""
@inline function get_secondary_code(galileo_e1b::GalileoE1B)
    1
end

"""
$(SIGNATURES)

Get the center (carrier) frequency for Galileo E1B.

# Returns
- `Frequency`: 1575.42 MHz

# Examples
```julia-repl
julia> get_center_frequency(GalileoE1B())
1575420000 Hz
```
"""
@inline function get_center_frequency(galileo_e1b::GalileoE1B)
    1_575_420_000Hz
end

"""
$(SIGNATURES)

Get the code chipping rate for Galileo E1B.

# Returns
- `Frequency`: 1.023 MHz

# Examples
```julia-repl
julia> get_code_frequency(GalileoE1B())
1023000 Hz
```
"""
@inline function get_code_frequency(galileo_e1b::GalileoE1B)
    1023_000Hz
end

"""
$(SIGNATURES)

Get the data bit rate for Galileo E1B.

# Returns
- `Frequency`: 250 Hz

# Examples
```julia-repl
julia> get_data_frequency(GalileoE1B())
250 Hz
```
"""
@inline function get_data_frequency(galileo_e1b::GalileoE1B)
    250Hz
end
