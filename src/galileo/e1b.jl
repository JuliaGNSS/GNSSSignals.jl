"""
    GalileoE1B{C} <: AbstractGalileoSignal{C}

Galileo E1B signal (the data-carrying component of Galileo E1 OS).

CBOC(6,1,1/11) modulation with a 4092-chip code at 1.023 Mcps, transmitted on
the E1 band — which shares the L1 carrier at 1575.42 MHz, so [`get_band`](@ref)
returns [`L1`](@ref) here.

# Example
```julia
e1b = GalileoE1B()
get_code_length(e1b)   # 4092
get_band(e1b)          # L1()
```
"""
struct GalileoE1B{C<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

get_modulation(::Type{<:GalileoE1B}) = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11)
@inline get_modulation(::GalileoE1B) = CBOC(BOCsin(1, 1), BOCsin(6, 1), 10 / 11)

"""
$(SIGNATURES)

Get the band the signal is transmitted on.

Galileo E1 shares the L1 carrier frequency (1575.42 MHz), so this returns
[`L1`](@ref) — band identity is by RF, not by ICD label.

# Examples
```julia-repl
julia> get_band(GalileoE1B())
L1()
```
"""
@inline get_band(::Type{<:GalileoE1B}) = L1()

"""
$(SIGNATURES)

Get the human-readable signal name.

# Examples
```julia-repl
julia> get_signal_name(GalileoE1B())
"Galileo E1B"
```
"""
get_signal_name(::GalileoE1B) = "Galileo E1B"

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
    codes = widen_codes_to_storage(read_galileo_e1b_codes())
    lut = build_signal_lut(get_modulation(GalileoE1B), codes, NoSecondaryCode())
    GalileoE1B(codes, lut)
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

Get the code chipping rate for Galileo E1B.

# Returns
- `Frequency`: 1.023 MHz

# Examples
```julia-repl
julia> get_code_frequency(GalileoE1B())
1023000 Hz
```
"""
@inline function get_code_frequency(::Type{<:GalileoE1B})
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
@inline function get_data_frequency(::Type{<:GalileoE1B})
    250Hz
end

"""
    GalileoE1B_BOC11{C} <: AbstractGalileoSignal{C}

BOC(1,1) approximation of Galileo E1B (the data-carrying component of
Galileo E1 OS).

Galileo E1B is specified as CBOC(6,1,1/11) — a Float32 weighted sum of
BOC(1,1) and BOC(6,1) requiring fs ≥ 2 · 6 · 1.023 MHz to fully
sample. Many software receivers substitute a pure BOC(1,1) replica
because (a) the BOC(6,1) component carries only 1/11 of the signal
power, so the correlation loss is ≈ 0.45 dB, and (b) BOC(1,1) needs
only fs ≥ 2 · 1.023 MHz, allowing lower sampling rates.
PocketSDR, for example, uses this substitution by default
(see `mod_code` in `sdr_code.c`).

Use this type when you want the lower sampling-rate variant; use
[`GalileoE1B`](@ref) for the full CBOC spec approximation. Either way
`gen_code!` emits `Int8` — only the single-chip [`get_code`](@ref)
accessor returns `Float32` (`get_code_type` is `Int16` here).

Primary code, code length, code frequency, data rate, and band are
identical to [`GalileoE1B`](@ref); only `get_modulation` differs.

# Example
```julia
e1b = GalileoE1B_BOC11()
get_modulation(e1b)    # BOCsin(1, 1)
get_code_length(e1b)   # 4092
```
"""
struct GalileoE1B_BOC11{C<:AbstractMatrix} <: AbstractGalileoSignal{C}
    codes::C
    lut::SignalLUT    # embedded per-signal LUT, always populated; see `build_signal_lut` / `gen_code!`
end

get_modulation(::Type{<:GalileoE1B_BOC11}) = BOCsin(1, 1)
@inline get_modulation(::GalileoE1B_BOC11) = BOCsin(1, 1)

@inline get_band(::Type{<:GalileoE1B_BOC11}) = L1()
get_signal_name(::GalileoE1B_BOC11) = "Galileo E1B (BOC(1,1) approximation)"

function GalileoE1B_BOC11()
    codes = widen_codes_to_storage(read_galileo_e1b_codes())
    lut = build_signal_lut(get_modulation(GalileoE1B_BOC11), codes, NoSecondaryCode())
    GalileoE1B_BOC11(codes, lut)
end

@inline get_code_length(::Type{<:GalileoE1B_BOC11}) = 4092
@inline get_secondary_code(::GalileoE1B_BOC11) = NoSecondaryCode()
@inline get_code_frequency(::Type{<:GalileoE1B_BOC11}) = 1023_000Hz
@inline get_data_frequency(::Type{<:GalileoE1B_BOC11}) = 250Hz
