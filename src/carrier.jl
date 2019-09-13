const CIS_LUT = SVector{64}(cis.((0:63) / 64 * 2π))
function cis_fast(x)
    @inbounds CIS_LUT[(floor(Int, x / 2π * 64) & 63) + 1]
end
