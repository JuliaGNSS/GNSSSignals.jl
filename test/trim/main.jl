# Build and run with the commands in the CI trim job. Safe trimming rejects
# unresolved calls even in error paths that valid signal construction never takes.
using GNSSSignals

function check_construction(::Type{S}, expected_factor::Int) where {S<:AbstractGNSSSignal}
    signal = S()
    valid = size(get_codes(signal), 1) == get_code_length(signal) &&
            signal.lut.subchip_factor == expected_factor &&
            signal.lut.table_length >= get_code_length(signal) * expected_factor &&
            isfinite(get_code(signal, 0.25, 1))
    if !valid
        println(Core.stderr, "constructed replica has inconsistent code metadata")
        return false
    end
    return true
end

function (@main)(args::Vector{String})::Cint
    valid = check_construction(GPSL1CA, 1)
    valid &= check_construction(GalileoE1B, 2)
    valid &= check_construction(GalileoE1B_CBOC, 12)
    valid &= check_construction(GalileoE1C, 2)
    valid &= check_construction(GalileoE1C_CBOC, 12)
    valid &= check_construction(GPSL1C_P, 2)
    valid &= check_construction(GPSL1C_P_TMBOC, 12)
    if valid
        println(Core.stdout, "default and full replica construction passed")
    end
    return valid ? 0 : 1
end
