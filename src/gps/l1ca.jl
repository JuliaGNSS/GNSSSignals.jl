"""
    GPSL1CA{C} <: AbstractGNSSSignal{C}

GPS L1 C/A signal.

The legacy civil GPS signal on the L1 band (1575.42 MHz). BPSK-modulated
1023-chip C/A code at 1.023 Mcps, fully modulated by the 50 bps navigation
message (no pilot component).

# Example
```julia
gpsl1ca = GPSL1CA()
get_code_length(gpsl1ca)   # 1023
get_band(gpsl1ca)          # L1()
```
"""
struct GPSL1CA{C<:AbstractMatrix} <: AbstractGNSSSignal{C}
    codes::C
end

get_modulation(::Type{<:GPSL1CA}) = LOC()
@inline get_modulation(::GPSL1CA) = LOC()

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

# Examples
```julia-repl
julia> get_band(GPSL1CA())
L1()
```
"""
@inline get_band(::GPSL1CA) = L1()

"""
$(SIGNATURES)

Get the human-readable signal name.

# Examples
```julia-repl
julia> get_signal_name(GPSL1CA())
"GPS L1 C/A"
```
"""
get_signal_name(::GPSL1CA) = "GPS L1 C/A"

function read_gpsl1ca_codes()
    read_in_codes(
        Int8,
        joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_gps_l1ca.bin"),
        37,
        1023,
    )
end

function GPSL1CA()
    GPSL1CA(widen_codes_to_storage(read_gpsl1ca_codes()))
end

"""
$(SIGNATURES)

Get the code length for GPS L1 C/A.

# Returns
- `Int`: 1023 chips

# Examples
```julia-repl
julia> get_code_length(GPSL1CA())
1023
```
"""
@inline function get_code_length(::GPSL1CA)
    1023
end

"""
$(SIGNATURES)

Get the secondary code for GPS L1 C/A.

GPS L1 C/A has no secondary code.

# Returns
- [`NoSecondaryCode`](@ref)

# Examples
```julia-repl
julia> get_secondary_code(GPSL1CA())
NoSecondaryCode()
```
"""
@inline function get_secondary_code(::GPSL1CA)
    NoSecondaryCode()
end

"""
$(SIGNATURES)

Get the code chipping rate for GPS L1 C/A.

# Returns
- `Frequency`: 1.023 MHz

# Examples
```julia-repl
julia> get_code_frequency(GPSL1CA())
1023000 Hz
```
"""
@inline function get_code_frequency(::GPSL1CA)
    1_023_000Hz
end

"""
$(SIGNATURES)

Get the data bit rate for GPS L1 C/A.

# Returns
- `Frequency`: 50 Hz

# Examples
```julia-repl
julia> get_data_frequency(GPSL1CA())
50 Hz
```
"""
@inline function get_data_frequency(::GPSL1CA)
    50Hz
end
