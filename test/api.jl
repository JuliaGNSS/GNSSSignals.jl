# Keeps the API reference, across its three pages, honest: every export is documented,
# every documented name exists, and the documented names that are not exported are
# exactly the ones meant to be reached qualified.
@testset "API reference" begin
    documented = Set{Symbol}()
    for page in ("api.md", "bands_and_time_systems.md", "codes.md"),
        line in eachline(joinpath(@__DIR__, "..", "docs", "src", page))

        m = match(r"^GNSSSignals\.(\w+!?)$", strip(line))
        isnothing(m) || push!(documented, Symbol(m[1]))
    end
    exported = Set(filter(!=(:GNSSSignals), names(GNSSSignals)))
    @test isempty(setdiff(exported, documented))
    @test all(name -> isdefined(GNSSSignals, name), documented)
    @test setdiff(documented, exported) == Set([
        :Modulation,
        :LOC,
        :BOC,
        :BOCsin,
        :BOCcos,
        :CBOC,
        :CodeFillEngine,
        :CodeFillState,
        :get_code_spectrum_BPSK,
        :get_code_spectrum_BOCsin,
        :get_code_spectrum_BOCcos,
        :read_in_codes,
        :SECONDS_PER_WEEK,
    ])
end
