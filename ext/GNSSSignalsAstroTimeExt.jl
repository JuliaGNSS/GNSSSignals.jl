module GNSSSignalsAstroTimeExt

using GNSSSignals: GNSSSignals, TAITime
using AstroTime: AstroTime, TAIEpoch

# AstroTime carries an additional compensation term; round it into our fraction.
AstroTime.TAIEpoch(t::TAITime) = TAIEpoch(t.second, t.fraction)
GNSSSignals.TAITime(epoch::TAIEpoch) = TAITime(epoch.second, epoch.fraction + epoch.error)

end
