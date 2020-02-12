const GPS_CA_CODES = extend_front_and_back(read_in_codes(
    joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_gps_l1.bin"),
    32,
    1023
))

"""
$(SIGNATURES)

Get codes of type GPSL1 as a Matrix where each column
represents a PRN.
```julia-repl
julia> get_code(GPSL1)
```
"""
function get_codes(::Type{GPSL1})
    @view GPS_CA_CODES[2:end - 2, :]
end

"""
$(SIGNATURES)

Get code length of GNSS system GPSL1.
```julia-repl
julia> get_code_length(GPSL1)
```
"""
@inline function get_code_length(::Type{GPSL1})
    1023
end

"""
$(SIGNATURES)

Get secondary code length of GNSS system GPSL1.
```julia-repl
julia> get_secondary_code_length(GPSL1)
```
"""
@inline function get_secondary_code_length(::Type{GPSL1})
    1
end

"""
$(SIGNATURES)

Get center frequency of GNSS system GPSL1.
```julia-repl
julia> get_center_frequency(GPSL1)
```
"""
@inline function get_center_frequency(::Type{GPSL1})
    1_575_420_000Hz
end

"""
$(SIGNATURES)

Get code frequency of GNSS system GPSL1.
```julia-repl
julia> get_code_frequency(GPSL1)
```
"""
@inline function get_code_frequency(::Type{GPSL1})
    1_023_000Hz
end

"""
$(SIGNATURES)

Get data frequency of GNSS system GPSL1.
```julia-repl
julia> get_data_frequency(GPSL1)
```
"""
@inline function get_data_frequency(::Type{GPSL1})
    50Hz
end

"""
$(SIGNATURES)

Get code of GNSS system GPSL1 at phase `phase` of prn `prn`.
The phase will not be wrapped by the code length. The phase has to smaller
than the code length and must be an integer.
```julia-repl
julia> get_code_unsafe(GPSL1, 10, 1)
```
"""
Base.@propagate_inbounds function get_code_unsafe(::Type{GPSL1}, phase::Int, prn::Int)
    GPS_CA_CODES[2 + phase, prn]
end
