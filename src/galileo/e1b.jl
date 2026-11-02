"""
    GalileoE1B{C, R} <: AbstractGalileoSignal{C}
    GalileoE1B()
    GalileoE1B(BOCsin)
    GalileoE1B(CBOC)

Galileo E1B signal (the data-carrying component of Galileo E1 OS).

Broadcast as CBOC(6,1,1/11) with a 4092-chip code at 1.023 Mcps, on the E1 band —
which shares the L1 carrier at 1575.42 MHz, so [`get_band`](@ref) returns [`L1`](@ref)
here.

The replica's modulation is chosen at construction and carried in the type parameter `R`:

- `GalileoE1B()`, the same as `GalileoE1B(BOCsin)`, replicates the BOC(1,1) component
  only. CBOC(6,1,1/11) is a weighted sum of BOC(1,1) and BOC(6,1) that needs
  fs ≥ 2 · 6 · 1.023 MHz to sample fully. Many software receivers substitute a pure
  BOC(1,1) replica because (a) the BOC(6,1) component carries only 1/11 of the signal
  power, so the correlation loss is ≈ 0.41 dB (`10 log10(10/11)`) against an ideal
  wideband signal, and (b) BOC(1,1) needs only fs ≥ 2 · 1.023 MHz. PocketSDR's
  [Python E1B generator](https://github.com/tomojitakasu/PocketSDR/blob/eda83016229dbc1374f8db4be81dc52f0f07bc13/python/sdr_code.py#L1436-L1443)
  uses this approximation; its
  [C E1B generator](https://github.com/tomojitakasu/PocketSDR/blob/eda83016229dbc1374f8db4be81dc52f0f07bc13/src/sdr_code.c#L1408-L1421)
  uses full CBOC.
- `GalileoE1B(CBOC)` replicates the full CBOC(6,1,1/11) signal.

Both are replicas of one transmitted component: they share the primary code, code length,
chip rate, bit rate, band, relative power, signal id (`:GalileoE1B`) and name.
[`get_modulation`](@ref) tells them apart. [`gen_code!`](@ref) emits `Int8` for both;
the single-chip [`get_code`](@ref) returns `Int` for BOC(1,1) and `Float32` for CBOC.

# Example
```julia
e1b = GalileoE1B()
get_code_length(e1b)                  # 4092
get_band(e1b)                         # L1()
get_modulation(e1b)                   # BOCsin(1, 1)
get_modulation(GalileoE1B(CBOC))      # CBOC(BOCsin(1, 1), BOCsin(6, 1), 10/11)
```
"""
struct GalileoE1B{C<:AbstractMatrix,R<:Modulation} <: AbstractGalileoSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

get_modulation(::Type{<:GalileoE1B{<:AbstractMatrix,BOCsin}}) = BOCsin(1, 1)
# The bare type stands for the default replica, as `GalileoE1B()` does.
get_modulation(::Type{GalileoE1B}) = BOCsin(1, 1)
get_modulation(::Type{<:GalileoE1B{<:AbstractMatrix,CBOC}}) =
    CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11)

@inline get_band(::Type{<:GalileoE1B}) = L1()

# 50/50 E1-B/E1-C power sharing of the E1 composite, whose −157.25 dBW total
# (Galileo OS SIS ICD v2.2, Table 13) is the unit here. This is the split *between*
# the two components; the 10/11 CBOC split *within* each one lives in
# [`get_modulation`](@ref). See [`get_relative_power`](@ref). A BOC(1,1) replica stands
# in for the same transmitted component and so carries its power: dropping the BOC(6,1)
# content changes the correlation shape, not the split between E1-B and E1-C.
@inline get_relative_power(::Type{<:GalileoE1B}) = 0.5

@inline get_signal_name(::Type{<:GalileoE1B}) = "Galileo E1B"

function read_galileo_e1b_codes()
    read_in_codes(
        Int8,
        joinpath(dirname(pathof(GNSSSignals)), "..", "data", "codes_galileo_e1b.bin"),
        50,
        4092,
    )
end

GalileoE1B() = GalileoE1B(BOCsin)
GalileoE1B(::Type{BOCsin}) = _galileo_e1b(BOCsin)
GalileoE1B(::Type{CBOC}) = _galileo_e1b(CBOC)

function _galileo_e1b(::Type{R}) where {R<:Modulation}
    codes = read_galileo_e1b_codes()
    S = GalileoE1B{typeof(codes),R}
    S(codes, build_signal_lut(get_modulation(S), codes, NoSecondaryCode()))
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
@inline function get_code_length(::Type{<:GalileoE1B})
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
julia> get_secondary_code(GalileoE1B())
NoSecondaryCode()
```
"""
@inline function get_secondary_code(::GalileoE1B)
    NoSecondaryCode()
end

"""
$(SIGNATURES)

Get the nominal chip rate of Galileo E1B.

# Returns
- `Frequency`: 1.023 MHz

# Examples
```julia-repl
julia> get_chip_rate(GalileoE1B())
1023000 Hz
```
"""
@inline function get_chip_rate(::Type{<:GalileoE1B})
    1023_000Hz
end

"""
$(SIGNATURES)

Get the bit rate, as tracked, of Galileo E1B.

The I/NAV message is tracked at 250 bit/s: the symbol rate of its rate-½ convolutional
code, the message carrying 125 bit/s.

# Returns
- `Frequency`: 250 Hz

# Examples
```julia-repl
julia> get_bit_rate(GalileoE1B())
250 Hz
```
"""
@inline function get_bit_rate(::Type{<:GalileoE1B})
    250Hz
end
