module GNSSSignalsAstroTimeExt

using GNSSSignals: GNSSSignals, TAITime
using AstroTime: AstroTime, TAIEpoch

# Both count whole TAI seconds since J2000 plus a fraction, so the conversions are exact.
AstroTime.TAIEpoch(t::TAITime) = TAIEpoch(t.second, t.fraction)
GNSSSignals.TAITime(epoch::TAIEpoch) = TAITime(epoch.second, epoch.fraction)

end
