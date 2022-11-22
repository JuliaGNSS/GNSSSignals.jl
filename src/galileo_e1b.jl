"""
GalileoE1B
"""
struct GalileoE1B{C <: AbstractMatrix} <: AbstractGNSS{C}
    codes::C
end

get_modulation(::Type{<:GalileoE1B}) = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11)

get_system_string(s::GalileoE1B) = "GalileoE1B"

function read_from_documentation(raw_code)
    raw_code_without_spaces = replace(replace(raw_code, " " => ""), "\n" => "")
    code_hex_array = map(x -> parse(UInt16, x, base = 16), collect(raw_code_without_spaces))
    code_bit_string = string(string.(code_hex_array, base = 2, pad = 4)...)
    map(x -> parse(Int8, x, base = 2), collect(code_bit_string)) .* Int8(2) .- Int8(1)
end

function read_galileo_e1b_codes()
    read_in_codes(
        Int8,
        joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_galileo_e1b.bin"),
        50,
        4092
    )
end

function GalileoE1B(;use_gpu = Val(false))
    _GalileoE1B(use_gpu)
end

function _GalileoE1B(::Val{false})
    GalileoE1B(Int16.(read_galileo_e1b_codes()))
end
function _GalileoE1B(::Val{true})
    GalileoE1B(CuMatrix{Float32}(read_galileo_e1b_codes()))
end

"""
$(SIGNATURES)

Get code length of GNSS system GalileoE1B.
```julia-repl
julia> get_code_length(GalileoE1B)
```
"""
@inline function get_code_length(galileo_e1b::GalileoE1B)
    4092
end

"""
$(SIGNATURES)

Get shortest code length of GNSS system GalileoE1B.
```julia-repl
julia> get_shortest_code_length(GalileoE1B)
```
"""
@inline function get_secondary_code_length(galileo_e1b::GalileoE1B)
    1
end

"""
$(SIGNATURES)

Get center frequency of GNSS system GalileoE1B.
```julia-repl
julia> get_center_frequency(GalileoE1B)
```
"""
@inline function get_center_frequency(galileo_e1b::GalileoE1B)
    1_575_420_000Hz
end

"""
$(SIGNATURES)

Get code frequency of GNSS system GalileoE1B.
```julia-repl
julia> get_code_frequency(GalileoE1B)
```
"""
@inline function get_code_frequency(galileo_e1b::GalileoE1B)
    1023_000Hz
end

"""
$(SIGNATURES)

Get data frequency of GNSS system GalileoE1B.
```julia-repl
julia> get_data_frequency(GalileoE1B)
```
"""
@inline function get_data_frequency(galileo_e1b::GalileoE1B)
    250Hz
end
