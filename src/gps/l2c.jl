"""
    GPSL2CM{C} <: AbstractGPSSignal{C}

GPS L2 CM signal (the moderate-length, data-carrying component of the GPS L2
civil signal L2C).

BPSK-modulated 10230-chip code at 511.5 kcps on the L2 band (1227.6 MHz),
giving a 20 ms code period (IS-GPS-200N §3.2.1.4). It carries the CNAV
message: a 25 bps data stream, rate-½ convolutionally encoded to a 50 sps
symbol stream (IS-GPS-200N §3.2.2); `get_data_frequency` returns the broadcast
symbol rate, matching the convention used across this package (see
[`GPSL5I`](@ref), [`GPSL1C_D`](@ref)).

In the broadcast L2C signal the L2 CM-code is time-division multiplexed
chip-by-chip with the dataless [`GPSL2CL`](@ref) pilot at 1.023 MHz; like
GNSS-SDR and PocketSDR this implementation models the CM component on its own
at its native 511.5 kcps chip rate, which is what a receiver replicates when
acquiring and tracking L2 CM. There is no secondary/overlay code.

The CM and CL codes share one degree-27 modular shift-register generator
(IS-GPS-200N §3.3.2.4, `gen_l2c_code`) and differ only in the per-PRN initial
register state (`GPS_L2CM_INITIAL_STATES`) and the short-cycle period. PRNs
1-63 are supported.

# Example
```julia
gpsl2cm = GPSL2CM()
get_code_length(gpsl2cm)      # 10230
get_code_frequency(gpsl2cm)   # 511500 Hz
get_band(gpsl2cm)             # L2()
```
"""
struct GPSL2CM{C<:AbstractMatrix} <: AbstractGPSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
    GPSL2CL{C} <: AbstractGPSSignal{C}

GPS L2 CL signal (the long, dataless pilot component of the GPS L2 civil
signal L2C).

BPSK-modulated 767250-chip code at 511.5 kcps on the L2 band (1227.6 MHz),
giving a 1.5 s code period (IS-GPS-200N §3.2.1.5). As the pilot component it
carries no navigation data, so [`get_data_frequency`](@ref) returns 0 Hz, and
it has no secondary/overlay code (the 1.5 s code is itself the long code).

In the broadcast L2C signal the dataless CL-code is time-division multiplexed
chip-by-chip with the [`GPSL2CM`](@ref) data component at 1.023 MHz; like
GNSS-SDR and PocketSDR this implementation models the CL component on its own
at its native 511.5 kcps chip rate.

The CL code uses the same degree-27 modular shift-register generator as
[`GPSL2CM`](@ref) (IS-GPS-200N §3.3.2.4, `gen_l2c_code`) and differs only in
the per-PRN initial register state (`GPS_L2CL_INITIAL_STATES`) and the longer
767250-chip short-cycle period. PRNs 1-63 are supported.

# Example
```julia
gpsl2cl = GPSL2CL()
get_code_length(gpsl2cl)      # 767250
get_data_frequency(gpsl2cl)   # 0 Hz
get_band(gpsl2cl)             # L2()
```
"""
struct GPSL2CL{C<:AbstractMatrix} <: AbstractGPSSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

"""
$(SIGNATURES)

Generate one GPS L2 civil code (CM or CL) from its initial shift-register
state `initial_state`, returning a length-`code_length` `Vector{Int8}` of ±1
chips.

The L2 CM- and L2 CL-codes share a single degree-27 modular (Galois)
shift-register generator (IS-GPS-200N §3.3.2.4); they differ only in the
per-PRN initial state and the number of chips before the register is
short-cycled (10230 for CM, 767250 for CL). At each step the register outputs
its least-significant bit (0 → +1, 1 → −1, the IS-GPS-200N / PocketSDR chip
convention), then shifts right and XORs the feedback mask
`GPS_L2C_FEEDBACK_MASK` back in when the output bit was 1.

```julia-repl
julia> code = gen_l2c_code(GPS_L2CM_INITIAL_STATES[1], GPS_L2CM_CODE_LENGTH);

julia> length(code)
10230
```
"""
function gen_l2c_code(initial_state, code_length)
    register = UInt32(initial_state)
    code = Vector{Int8}(undef, code_length)
    @inbounds for i = 1:code_length
        output = register & one(register)
        code[i] = iszero(output) ? Int8(1) : Int8(-1)
        register = (register >> 1) ⊻ (GPS_L2C_FEEDBACK_MASK * output)
    end
    return code
end

# Build the `(code_length, GPS_L2C_NUM_PRNS)` primary-code matrix from a
# per-PRN initial-state table, preallocated and filled column-by-column (the
# pattern used by `_l1c_build_primary_codes`).
function _build_l2c_codes(initial_states::AbstractVector, code_length::Integer)
    codes = Matrix{Int8}(undef, code_length, GPS_L2C_NUM_PRNS)
    for prn = 1:GPS_L2C_NUM_PRNS
        codes[:, prn] = gen_l2c_code(initial_states[prn], code_length)
    end
    return codes
end

read_gpsl2cm_codes() = _build_l2c_codes(GPS_L2CM_INITIAL_STATES, GPS_L2CM_CODE_LENGTH)
read_gpsl2cl_codes() = _build_l2c_codes(GPS_L2CL_INITIAL_STATES, GPS_L2CL_CODE_LENGTH)

function GPSL2CM()
    codes = widen_codes_to_storage(read_gpsl2cm_codes())
    lut = build_signal_lut(get_modulation(GPSL2CM), codes, NoSecondaryCode())
    GPSL2CM(codes, lut)
end

function GPSL2CL()
    codes = widen_codes_to_storage(read_gpsl2cl_codes())
    lut = build_signal_lut(get_modulation(GPSL2CL), codes, NoSecondaryCode())
    GPSL2CL(codes, lut)
end

# Shared interface (modulation, band, frequencies).

get_modulation(::Type{<:GPSL2CM}) = LOC()
@inline get_modulation(::GPSL2CM) = LOC()
get_modulation(::Type{<:GPSL2CL}) = LOC()
@inline get_modulation(::GPSL2CL) = LOC()

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

# Examples
```julia-repl
julia> get_band(GPSL2CM())
L2()
```
"""
@inline get_band(::Type{<:GPSL2CM}) = L2()
@inline get_band(::Type{<:GPSL2CL}) = L2()

# L2C total −160.0 dBW (IS-GPS-200N, Table 3-Va, IIR-M/IIF worst case; the ICD
# gives the combined L2C only). CM and CL are chip-by-chip time-multiplexed, so
# each carries half the average power: −160.0 + 10·log10(0.5) = −163.0 dBW each.
# This split is derived, not tabulated. See [`get_min_received_power`](@ref).
@inline get_min_received_power(::Type{<:GPSL2CM}) = _dbw_to_watts(-160.0 + 10log10(0.5))
@inline get_min_received_power(::Type{<:GPSL2CL}) = _dbw_to_watts(-160.0 + 10log10(0.5))

"""
$(SIGNATURES)

Get the human-readable signal name.

# Examples
```julia-repl
julia> get_signal_name(GPSL2CM())
"GPS L2CM"
```
"""
get_signal_name(::GPSL2CM) = "GPS L2CM"
get_signal_name(::GPSL2CL) = "GPS L2CL"

"""
$(SIGNATURES)

Get the code length for GPS L2 CM (10230 chips, 20 ms at 511.5 kcps).

# Examples
```julia-repl
julia> get_code_length(GPSL2CM())
10230
```
"""
@inline get_code_length(::Type{<:GPSL2CM}) = GPS_L2CM_CODE_LENGTH

"""
$(SIGNATURES)

Get the code length for GPS L2 CL (767250 chips, 1.5 s at 511.5 kcps).

# Examples
```julia-repl
julia> get_code_length(GPSL2CL())
767250
```
"""
@inline get_code_length(::Type{<:GPSL2CL}) = GPS_L2CL_CODE_LENGTH

"""
$(SIGNATURES)

Get the secondary code for GPS L2 CM.

The L2 CM-code has no secondary/overlay code.

# Returns
- [`NoSecondaryCode`](@ref)
"""
@inline get_secondary_code(::GPSL2CM) = NoSecondaryCode()

"""
$(SIGNATURES)

Get the secondary code for GPS L2 CL.

The L2 CL-code has no secondary/overlay code; the 1.5 s code is itself the
long code.

# Returns
- [`NoSecondaryCode`](@ref)
"""
@inline get_secondary_code(::GPSL2CL) = NoSecondaryCode()

"""
$(SIGNATURES)

Get the code chipping rate for the GPS L2 civil signals (511.5 kHz, both
components).

# Examples
```julia-repl
julia> get_code_frequency(GPSL2CM())
511500 Hz
```
"""
@inline get_code_frequency(::Type{<:GPSL2CM}) = GPS_L2C_CODE_FREQUENCY
@inline get_code_frequency(::Type{<:GPSL2CL}) = GPS_L2C_CODE_FREQUENCY

"""
$(SIGNATURES)

Get the data symbol rate for GPS L2 CM.

The L2 CM-code carries the CNAV message at 50 sps (25 bps with rate-½
convolutional coding, IS-GPS-200N §3.2.2); `get_data_frequency` returns the
broadcast symbol rate, matching the convention used across this package (see
[`GPSL5I`](@ref)).

# Returns
- `Frequency`: 50 Hz
"""
@inline get_data_frequency(::Type{<:GPSL2CM}) = 50Hz

"""
$(SIGNATURES)

Get the data symbol rate for GPS L2 CL.

The L2 CL-code is a dataless pilot, so its data frequency is 0 Hz.

# Returns
- `Frequency`: 0 Hz
"""
@inline get_data_frequency(::Type{<:GPSL2CL}) = 0Hz
