using AeroPanels
using Test

@testset "AeroPanels.jl" begin
    include("drag/DragVerification.jl")
    include("sweep/SweepVerification.jl")
    include("fmu/test_fmu_construction.jl")
end
