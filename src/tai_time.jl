# A minimal, trim-safe TAI epoch, standing in for AstroTime.jl's `TAIEpoch`, whose
# `__init__` a `juliac --trim=safe` build cannot verify (it registers `Dates` format tokens
# and, through EarthOrientation and RemoteFiles, looks for data files and download tools).
# With AstroTime loaded, the `GNSSSignalsAstroTimeExt` extension converts between the two.

# J2000, as a TAI calendar label: the origin `second` counts from, as in AstroTime.
const J2000 = DateTime(2000, 1, 1, 12)

"""
    TAITime(second::Integer, fraction::Real)
    TAITime(datetime::Dates.DateTime)

An epoch on International Atomic Time (TAI): `second` whole seconds after J2000
(`2000-01-01T12:00:00` TAI) plus `fraction ∈ [0, 1)` of a second, AstroTime's `TAIEpoch`
layout, so the two convert exactly (with AstroTime loaded, `TAIEpoch(t)` and
`TAITime(epoch)`). A `fraction` outside `[0, 1)` is normalised into it. The time systems
date their weeks in it ([`get_week_start_epoch`](@ref)).

`TAITime(datetime)` reads a `DateTime` as a TAI calendar label (millisecond resolution),
such as [`get_tai_system_start_time`](@ref) returns; `Dates.DateTime(t)` gives the TAI
calendar label back, rounded to the millisecond. The difference of two `TAITime`s is in
seconds, `t + Δ` (or `Δ + t`) shifts `t` by `Δ` seconds, and `a ≈ b` holds when the
two are within `atol` seconds of each other (default `1e-9`, a third of a metre of range).

# Fields
- `second::Int64`: whole TAI seconds since J2000.
- `fraction::Float64`: the fraction of the second, in `[0, 1)`.

```julia-repl
julia> TAITime(get_tai_system_start_time(GPST()))
1980-01-06T00:00:19 TAI

julia> TAITime(0, 0.25) + 1.5
2000-01-01T12:00:01.750 TAI
```
"""
struct TAITime
    second::Int64
    fraction::Float64
    function TAITime(second::Integer, fraction::Real)
        whole = floor(Int64, fraction)
        # `+ 0.0` turns a `-0.0` fraction into `0.0`, so that equal epochs compare equal.
        rest = Float64(fraction) - whole + 0.0
        # A fraction just below a whole number rounds up to exactly `1.0` in the
        # subtraction (`-1e-17 - (-1)` is `1.0` in `Float64`), which would escape the
        # `[0, 1)` invariant; carry it into the whole seconds instead.
        rest == 1.0 && (whole += 1; rest = 0.0)
        new(Int64(second) + whole, rest)
    end
end

function TAITime(datetime::DateTime)
    milliseconds = Dates.value(datetime - J2000)
    TAITime(fld(milliseconds, 1000), mod(milliseconds, 1000) / 1000)
end

Dates.DateTime(t::TAITime) =
    J2000 + Dates.Millisecond(1000 * t.second + round(Int64, 1000 * t.fraction))

Base.:-(a::TAITime, b::TAITime) = (a.second - b.second) + (a.fraction - b.fraction)

function Base.:+(t::TAITime, Δ::Real)
    whole = floor(Int64, Δ)
    TAITime(t.second + whole, t.fraction + (Δ - whole))
end
Base.:+(Δ::Real, t::TAITime) = t + Δ
Base.:-(t::TAITime, Δ::Real) = t + (-Δ)

Base.isless(a::TAITime, b::TAITime) = isless((a.second, a.fraction), (b.second, b.fraction))

# On the difference, with an absolute tolerance. AstroTime compares the seconds since J2000
# with the default relative tolerance instead, which at today's ~8e8 s takes epochs about
# ten seconds apart for equal.
Base.isapprox(a::TAITime, b::TAITime; atol::Real = 1e-9) = abs(a - b) <= atol

Base.show(io::IO, t::TAITime) = print(io, DateTime(t), " TAI")
