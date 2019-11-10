function read_from_documentation(raw_code)
    raw_code_without_spaces = replace(replace(raw_code, " " => ""), "\n" => "")
    code_hex_array = map(x -> parse(UInt16, x, base = 16), collect(raw_code_without_spaces))
    code_bit_string = string(string.(code_hex_array, base = 2, pad = 4)...)
    map(x -> parse(Int8, x, base = 2), collect(code_bit_string)) .* Int8(2) .- Int8(1)
end

const GALILEO_E1B_CODES = read_in_codes(
    joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_galileo_e1b.bin"),
    50,
    4092
)

"""
$(SIGNATURES)

Get codes of type GalileoE1B as a Matrix where each column
represents a PRN.
```julia-repl
julia> get_code(GalileoE1B)
```
"""
function get_codes(::Type{GalileoE1B})
    GALILEO_E1B_CODES
end

"""
$(SIGNATURES)

Get code length of GNSS system GalileoE1B.
```julia-repl
julia> get_code_length(GalileoE1B)
```
"""
@inline function get_code_length(::Type{GalileoE1B})
    4092
end

"""
$(SIGNATURES)

Get shortest code length of GNSS system GalileoE1B.
```julia-repl
julia> get_shortest_code_length(GalileoE1B)
```
"""
@inline function get_secondary_code_length(::Type{GalileoE1B})
    1
end

"""
$(SIGNATURES)

Get center frequency of GNSS system GalileoE1B.
```julia-repl
julia> get_center_frequency(GalileoE1B)
```
"""
@inline function get_center_frequency(::Type{GalileoE1B})
    1_575_420_000Hz
end

"""
$(SIGNATURES)

Get code frequency of GNSS system GalileoE1B.
```julia-repl
julia> get_code_frequency(GalileoE1B)
```
"""
@inline function get_code_frequency(::Type{GalileoE1B})
    1023_000Hz
end

"""
$(SIGNATURES)

Get data frequency of GNSS system GalileoE1B.
```julia-repl
julia> get_data_frequency(GalileoE1B)
```
"""
@inline function get_data_frequency(::Type{GalileoE1B})
    250Hz
end

"""
$(SIGNATURES)

Get code of type GalileoE1B at phase `phase` of PRN `prn`. Includes only BOC(1,1) at the
moment.
```julia-repl
julia> get_code(GalileoE1B, 10.3, 1)
```
"""
Base.@propagate_inbounds function get_code(::Type{GalileoE1B}, phase, prn::Int)
    floored_2phase = floor(Int, 2 * phase)
    get_code_unsafe(
        GalileoE1B,
        mod(
            floored_2phase >> 1,
            get_code_length(GalileoE1B) * get_secondary_code_length(GalileoE1B)
        ),
        prn
    ) * (iseven(floored_2phase) * 2 - 1)
end

"""
$(SIGNATURES)

Get code of type GalileoE1B at phase `phase` of PRN `prn`.
The phase will not be wrapped by the code length. The phase has to smaller
than the code length. Includes only BOC(1,1) at the moment.
```julia-repl
julia> get_code_unsafe(GalileoE1B, 10.3, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(::Type{GalileoE1B}, phase, prn::Int)
    floored_2phase = floor(Int, 2 * phase)
    get_code_unsafe(GalileoE1B, floored_2phase >> 1, prn) *
        (iseven(floored_2phase) * 2 - 1)
end

"""
$(SIGNATURES)

Get code of GNSS system GalileoE1B at phase `phase` of prn `prn`.
The phase will not be wrapped by the code length. The phase has to smaller
than the code length and must be an integer.
```julia-repl
julia> get_code_unsafe(GalileoE1B, 10, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(::Type{GalileoE1B}, phase::Int, prn::Int)
    GALILEO_E1B_CODES[1 + phase, prn]
end
