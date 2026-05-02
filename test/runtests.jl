using AeroPanels
using Test

@testset "AeroPanels.jl" begin
    include("drag/DragVerification.jl")
    include("sweep/SweepVerification.jl")
end
