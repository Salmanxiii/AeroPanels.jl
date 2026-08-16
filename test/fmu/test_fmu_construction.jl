using AeroPanels
using Test
using LinearAlgebra
using ForwardDiff

@testset "FMU Model & Array Construction" begin
    # 1. Instantiate test model
    model = SE2ASteady()

    # 2. Test Idiomatic CreateFMUArray Dispatches
    @testset "CreateFMUArray Dispatches" begin
        arr_default = CreateFMUArray(model)
        @test arr_default isa AbstractFMUArray{Float64}
        @test length(arr_default.data) == GetTotalSize(model)

        arr_f32 = CreateFMUArray(Float32, model)
        @test arr_f32 isa AbstractFMUArray{Float32}
        @test length(arr_f32.data) == GetTotalSize(model)

        arr_f64 = CreateFMUArray(Float64, model)
        @test arr_f64 isa AbstractFMUArray{Float64}
        @test length(arr_f64.data) == GetTotalSize(model)

        custom_data = zeros(Float64, GetTotalSize(model))
        arr_custom = CreateFMUArray(model, custom_data)
        @test arr_custom isa AbstractFMUArray{Float64}
        @test pointer(arr_custom.data) == pointer(custom_data)
    end

    # 3. Test Metadata & Layout Reflection
    @testset "FMU Layout Reflection" begin
        layout = GetFMULayout(model)
        @test layout isa FMULayout
        @test NumberOfStates(model) >= 0
        @test GetTotalSize(model) > 0
    end

    # 4. Test Low-Level Solvers & AD ForwardDiff Compatibility
    @testset "Low-Level Solvers & AD Support" begin
        arr = CreateFMUArray(model)
        cache = AllocateFMUCaches(model)
        
        # Test Initialize & Evaluation
        InitializeFMU!(arr, cache, model, 0.0)
        EvaluateDerivatives!(arr, cache, model, 0.0)
        EvaluateOutputs!(arr, cache, model, 0.0)
        
        # Test Dual Number Buffer (Automatic Differentiation)
        dual_data = [ForwardDiff.Dual(arr.data[i], 1.0) for i in 1:length(arr.data)]
        dual_arr = CreateFMUArray(model, dual_data)
        @test dual_arr isa AbstractFMUArray{<:ForwardDiff.Dual}
    end
end
