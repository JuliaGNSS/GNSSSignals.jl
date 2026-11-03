"""
$(SIGNATURES)

Get the spectral power density of a GNSS signal at a given frequency.

Computes the power spectral density from the signal's modulation type, at baseband
frequency `f`.

# Examples
```julia-repl
julia> using Unitful: kHz
julia> get_code_spectrum(GPSL1CA(), 0kHz)
9.775171065493646e-7
```
"""
get_code_spectrum(signal, f) = get_code_spectrum(get_modulation(signal), signal, f)
get_code_spectrum(modulation::BPSK, signal, f) =
    get_code_spectrum_BPSK(get_chip_rate(signal), f)
# The chip rate f_c is the signal's (`n · 1.023 MHz` already) and the subcarrier frequency
# f_s = m · 1.023 MHz is `m/n` of it.
get_code_spectrum(modulation::BOCsin, signal, f) = get_code_spectrum_BOCsin(
    get_chip_rate(signal),
    (modulation.m / modulation.n) * get_chip_rate(signal),
    f,
)
get_code_spectrum(modulation::BOCcos, signal, f) = get_code_spectrum_BOCcos(
    get_chip_rate(signal),
    (modulation.m / modulation.n) * get_chip_rate(signal),
    f,
)
function get_code_spectrum(modulation::CBOC, signal, f)
    get_code_spectrum(modulation.boc1, signal, f) * modulation.boc1_power +
    get_code_spectrum(modulation.boc2, signal, f) * (1 - modulation.boc1_power)
end

# For TMBOC, weight by the fraction of chips that use each BOC variant.
function get_code_spectrum(modulation::TMBOC{B1, B2, N}, signal, f) where {B1, B2, N}
    n_boc2 = count(modulation.pattern)
    boc2_frac = n_boc2 / N
    get_code_spectrum(modulation.boc1, signal, f) * (1 - boc2_frac) +
    get_code_spectrum(modulation.boc2, signal, f) * boc2_frac
end


"""
$(SIGNATURES)

Calculate the spectral power density of a BPSK modulated signal.

Evaluated at baseband frequency `f` for chip rate `chip_rate`.

# Examples
```julia-repl
julia> using Unitful: MHz, kHz
julia> get_code_spectrum_BPSK(1.023MHz, 0kHz)
9.775171065493646e-7
```
"""
function get_code_spectrum_BPSK(chip_rate::Frequency, f)
    return get_code_spectrum_BPSK(chip_rate / 1Hz, f)
end
function get_code_spectrum_BPSK(chip_rate, f::Frequency)
    return get_code_spectrum_BPSK(chip_rate, f / 1Hz)
end
function get_code_spectrum_BPSK(chip_rate::Frequency, f::Frequency)
    return get_code_spectrum_BPSK(chip_rate / 1Hz, f / 1Hz)
end
function get_code_spectrum_BPSK(chip_rate, f)
    return sinc(f / chip_rate)^2 / chip_rate
end

"""
$(SIGNATURES)

Calculate the spectral power density of a sine-phased BOC modulated signal.

Evaluated at baseband frequency `f` for chip rate `chip_rate` and subcarrier frequency
`subcarrier_frequency`.

A chip spans `k = 2 · subcarrier_frequency / chip_rate` subcarrier half-periods, which must
be a whole number (an `ArgumentError` otherwise). An even `k` puts whole subcarrier periods in
a chip and nulls the spectrum at `f = 0`; an odd `k` leaves half a period over, so it does not
(Betz, "Binary Offset Carrier Modulations for Radionavigation", NAVIGATION 48(4), 2001).
"""
function get_code_spectrum_BOCsin(chip_rate::Frequency, subcarrier_frequency::Frequency, f)
    return get_code_spectrum_BOCsin(chip_rate / 1Hz, subcarrier_frequency / 1Hz, f)
end
function get_code_spectrum_BOCsin(chip_rate, subcarrier_frequency, f::Frequency)
    return get_code_spectrum_BOCsin(chip_rate, subcarrier_frequency, f / 1Hz)
end
function get_code_spectrum_BOCsin(chip_rate::Frequency, subcarrier_frequency::Frequency, f::Frequency)
    return get_code_spectrum_BOCsin(chip_rate / 1Hz, subcarrier_frequency / 1Hz, f / 1Hz)
end
function get_code_spectrum_BOCsin(chip_rate, subcarrier_frequency, f)
    k = _half_periods_per_chip(chip_rate, subcarrier_frequency)
    if iseven(k)
        return ((sinc(f / chip_rate) * tan(pi * f / (2 * subcarrier_frequency)))^2 / chip_rate)
    else
        return (sinc(f / (2 * subcarrier_frequency)) * cospi(f / chip_rate) /
                (k * cospi(f / (2 * subcarrier_frequency))))^2 / chip_rate
    end
end

"""
$(SIGNATURES)

Calculate the spectral power density of a cosine-phased BOC modulated signal.

Evaluated at baseband frequency `f` for chip rate `chip_rate` and subcarrier frequency
`subcarrier_frequency`.

As for [`get_code_spectrum_BOCsin`](@ref GNSSSignals.get_code_spectrum_BOCsin), a chip must
span a whole number `k = 2 · subcarrier_frequency / chip_rate` of subcarrier half-periods,
and an even and an odd `k` take different forms.
"""
function get_code_spectrum_BOCcos(chip_rate::Frequency, subcarrier_frequency::Frequency, f)
    return get_code_spectrum_BOCcos(chip_rate / 1Hz, subcarrier_frequency / 1Hz, f)
end
function get_code_spectrum_BOCcos(chip_rate, subcarrier_frequency, f::Frequency)
    return get_code_spectrum_BOCcos(chip_rate, subcarrier_frequency, f / 1Hz)
end
function get_code_spectrum_BOCcos(chip_rate::Frequency, subcarrier_frequency::Frequency, f::Frequency)
    return get_code_spectrum_BOCcos(chip_rate / 1Hz, subcarrier_frequency / 1Hz, f / 1Hz)
end
function get_code_spectrum_BOCcos(chip_rate, subcarrier_frequency, f)
    k = _half_periods_per_chip(chip_rate, subcarrier_frequency)
    if iseven(k)
        return (2 * sinc(f / chip_rate) * sinpi(f / (4 * subcarrier_frequency))^2 /
                cospi(f / (2 * subcarrier_frequency)))^2 / chip_rate
    else
        return (sinpi(f / (4 * subcarrier_frequency)) * sinc(f / (4 * subcarrier_frequency)) *
                cospi(f / chip_rate) / (k * cospi(f / (2 * subcarrier_frequency))))^2 / chip_rate
    end
end

# The number of subcarrier half-periods `k = 2 f_s / f_c` a chip spans. The closed-form BOC
# spectra hold for a whole `k` only, and take one form for an even and another for an odd one.
function _half_periods_per_chip(chip_rate, subcarrier_frequency)
    k = 2 * subcarrier_frequency / chip_rate
    k_whole = round(Int, k)
    isapprox(k, k_whole) || throw(ArgumentError(
        "a BOC chip must span a whole number of subcarrier half-periods; got $k"))
    return k_whole
end
