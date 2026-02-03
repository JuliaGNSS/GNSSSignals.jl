"""
    GPSL1{C} <: AbstractGNSS{C}

GPS L1 C/A signal type.

GPS L1 uses BPSK (LOC) modulation with a 1023-chip C/A code at 1.023 Mcps,
transmitted on the L1 carrier frequency of 1575.42 MHz.

# Example
```julia
gpsl1 = GPSL1()
get_code_length(gpsl1)  # 1023
```
"""
struct GPSL1{C<:AbstractMatrix} <: AbstractGNSS{C}
    codes::C
end

get_modulation(::Type{<:GPSL1}) = LOC()
@inline get_modulation(::GPSL1) = LOC()

"""
$(SIGNATURES)

Get the system name as a string.

# Arguments
- `system`: A GNSS system instance

# Returns
- `String`: System identifier string

# Examples
```julia-repl
julia> get_system_string(GPSL1())
"GPSL1"
```
"""
get_system_string(s::GPSL1) = "GPSL1"

function read_gpsl1_codes()
    read_in_codes(
        Int8,
        joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_gps_l1.bin"),
        37,
        1023,
    )
end

function GPSL1()
    GPSL1(Int16.(read_gpsl1_codes()))
end

"""
$(SIGNATURES)

Get the code length for GPS L1.

# Returns
- `Int`: 1023 chips

# Examples
```julia-repl
julia> get_code_length(GPSL1())
1023
```
"""
@inline function get_code_length(gpsl1::GPSL1)
    1023
end

"""
$(SIGNATURES)

Get the secondary code for GPS L1.

GPS L1 has no secondary code, so this returns 1.

# Returns
- `Int`: 1 (no secondary code)

# Examples
```julia-repl
julia> get_secondary_code(GPSL1())
1
```
"""
@inline function get_secondary_code(gpsl1::GPSL1)
    1
end

"""
$(SIGNATURES)

Get the center (carrier) frequency for GPS L1.

# Returns
- `Frequency`: 1575.42 MHz

# Examples
```julia-repl
julia> get_center_frequency(GPSL1())
1575420000 Hz
```
"""
@inline function get_center_frequency(gpsl1::GPSL1)
    1_575_420_000Hz
end

"""
$(SIGNATURES)

Get the code chipping rate for GPS L1.

# Returns
- `Frequency`: 1.023 MHz

# Examples
```julia-repl
julia> get_code_frequency(GPSL1())
1023000 Hz
```
"""
@inline function get_code_frequency(gpsl1::GPSL1)
    1_023_000Hz
end

"""
$(SIGNATURES)

Get the data bit rate for GPS L1.

# Returns
- `Frequency`: 50 Hz

# Examples
```julia-repl
julia> get_data_frequency(GPSL1())
50 Hz
```
"""
@inline function get_data_frequency(gpsl1::GPSL1)
    50Hz
end
