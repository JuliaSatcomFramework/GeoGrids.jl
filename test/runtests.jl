using TestItemRunner

include("filtering.jl")
include("grids.jl")
include("interfaces.jl")
include("plots.jl")
include("tesselation.jl")
include("offsetting.jl")

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(GeoGrids; ambiguities=false)
    Aqua.test_ambiguities(GeoGrids)
end

@testitem "JET" tags=[:jet] begin
    using JET
    report_package("GeoGrids")
end

@run_package_tests verbose = true