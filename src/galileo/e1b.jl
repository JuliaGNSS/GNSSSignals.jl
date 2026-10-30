"""
    GalileoE1B_CBOC{C} <: AbstractGalileoSignal{C}

Galileo E1B data-carrying component with the full CBOC(6,1,1/11)
replica. The signal uses a 4092-chip primary code at 1.023 Mcps and no secondary code.
It is transmitted on the E1 band, which shares the L1 carrier.

The CBOC replica uses in-phase BOC(6,1) relative to BOC(1,1).
Use [`GalileoE1B`](@ref) for a BOC(1,1) approximation. Both share
the same primary code, secondary code, rates, and relative power. `get_code`
returns `Float32` for CBOC; `gen_code!` emits an `Int8` replica.
"""
struct GalileoE1B_CBOC{C<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

get_modulation(::Type{<:GalileoE1B_CBOC}) = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11)

@inline get_band(::Type{<:GalileoE1B_CBOC}) = L1()

# 50/50 E1-B/E1-C power sharing of the E1 composite, whose −157.25 dBW total
# (Galileo OS SIS ICD v2.2, Table 13) is the unit here. This is the split *between*
# the two components; the 10/11 CBOC split *within* each one lives in
# [`get_modulation`](@ref). See [`get_relative_power`](@ref).
@inline get_relative_power(::Type{<:GalileoE1B_CBOC}) = 0.5

@inline get_signal_name(::Type{<:GalileoE1B_CBOC}) = "Galileo E1B (CBOC)"

function read_galileo_e1b_codes()
    read_in_codes(
        Int8,
        joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_galileo_e1b.bin"),
        50,
        4092,
    )
end

function GalileoE1B_CBOC()
    codes = read_galileo_e1b_codes()
    lut = build_signal_lut(get_modulation(GalileoE1B_CBOC), codes, NoSecondaryCode())
    GalileoE1B_CBOC(codes, lut)
end

"""
$(SIGNATURES)

Get the code length for Galileo E1B.

# Returns
- `Int`: 4092 chips

# Examples
```julia-repl
julia> get_code_length(GalileoE1B_CBOC())
4092
```
"""
@inline function get_code_length(::Type{<:GalileoE1B_CBOC})
    4092
end

"""
$(SIGNATURES)

Get the secondary code for Galileo E1B.

Galileo E1B has no secondary code (it is the data component; the secondary
code lives on the E1C pilot).

# Returns
- [`NoSecondaryCode`](@ref)

# Examples
```julia-repl
julia> get_secondary_code(GalileoE1B_CBOC())
NoSecondaryCode()
```
"""
@inline function get_secondary_code(::GalileoE1B_CBOC)
    NoSecondaryCode()
end

"""
$(SIGNATURES)

Get the nominal chip rate of Galileo E1B.

# Returns
- `Frequency`: 1.023 MHz

# Examples
```julia-repl
julia> get_chip_rate(GalileoE1B_CBOC())
1023000 Hz
```
"""
@inline function get_chip_rate(::Type{<:GalileoE1B_CBOC})
    1023_000Hz
end

"""
$(SIGNATURES)

Get the bit rate, as tracked, of Galileo E1B.

The I/NAV message is transmitted at 250 symbols/s: the symbol rate of its rate-½ convolutional
code, the message carrying 125 bit/s.

# Returns
- `Frequency`: 250 Hz

# Examples
```julia-repl
julia> get_bit_rate(GalileoE1B_CBOC())
250 Hz
```
"""
@inline function get_bit_rate(::Type{<:GalileoE1B_CBOC})
    250Hz
end

"""
    GalileoE1B{C} <: AbstractGalileoSignal{C}

BOC(1,1) replica of the Galileo E1B data-carrying component. The broadcast
waveform uses CBOC(6,1,1/11); use [`GalileoE1B_CBOC`](@ref) for the
full CBOC replica. Both variants use the same 4092-chip primary code at
1.023 Mcps and no secondary code.

The BOC(1,1) approximation omits the BOC(6,1) component, which carries 1/11
of the broadcast power. `get_code` and `gen_code!` return `Int8` values.
"""
struct GalileoE1B{C<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

get_modulation(::Type{<:GalileoE1B}) = BOCsin(1, 1)

@inline get_band(::Type{<:GalileoE1B}) = L1()

# Approximates E1-B, so it stands in for the same transmitted component and
# carries its power. The approximation drops the BOC(6,1) content, which changes
# the correlation shape, not the E1-B/E1-C power split this reports.
@inline get_relative_power(::Type{<:GalileoE1B}) = get_relative_power(GalileoE1B_CBOC)

@inline get_signal_name(::Type{<:GalileoE1B}) = "Galileo E1B (BOC(1,1) approximation)"

function GalileoE1B()
    codes = read_galileo_e1b_codes()
    lut = build_signal_lut(get_modulation(GalileoE1B), codes, NoSecondaryCode())
    GalileoE1B(codes, lut)
end

@inline get_code_length(::Type{<:GalileoE1B}) = 4092
@inline get_secondary_code(::GalileoE1B) = NoSecondaryCode()
@inline get_chip_rate(::Type{<:GalileoE1B}) = 1023_000Hz
@inline get_bit_rate(::Type{<:GalileoE1B}) = 250Hz
