struct GPSL1{C <: AbstractMatrix} <: AbstractGNSS
    codes::C
end

function GPSL1()
    codes = read_in_codes(
        joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_gps_l1.bin"),
        37,
        1023
    )
    GPSL1(extend_front_and_back(codes, size(codes, 1)))
end

"""
$(SIGNATURES)

Get code length of GNSS system GPSL1.
```julia-repl
julia> get_code_length(gpsl1)
```
"""
@inline function get_code_length(gpsl1::GPSL1)
    1023
end

"""
$(SIGNATURES)

Get secondary code length of GNSS system GPSL1.
```julia-repl
julia> get_secondary_code_length(gpsl1)
```
"""
@inline function get_secondary_code_length(gpsl1::GPSL1)
    1
end

"""
$(SIGNATURES)

Get center frequency of GNSS system GPSL1.
```julia-repl
julia> get_center_frequency(gpsl1)
```
"""
@inline function get_center_frequency(gpsl1::GPSL1)
    1_575_420_000Hz
end

"""
$(SIGNATURES)

Get code frequency of GNSS system GPSL1.
```julia-repl
julia> get_code_frequency(gpsl1)
```
"""
@inline function get_code_frequency(gpsl1::GPSL1)
    1_023_000Hz
end

"""
$(SIGNATURES)

Get data frequency of GNSS system GPSL1.
```julia-repl
julia> get_data_frequency(gpsl1)
```
"""
@inline function get_data_frequency(gpsl1::GPSL1)
    50Hz
end