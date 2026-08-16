# src/Interfaces/Simulation.jl

using OrdinaryDiffEqTsit5
using DiffEqCallbacks
using StructArrays
using StaticArrays

function BuildOutputTuple(array::AbstractFMUArray, model::AbstractFMUModel, saveVars)
    layout = GetFMULayout(model)
    out_states = hasproperty(layout, :OutputStates) ? layout.OutputStates : ()
    out_derivs = hasproperty(layout, :OutputDerivatives) ? layout.OutputDerivatives : ()
    all_outputs = (out_states..., out_derivs..., layout.Outputs...)
    
    selected_vars = if saveVars == :summary
        filter(v -> v != :nodalForces, all_outputs)
    elseif saveVars == :all
        all_outputs
    elseif saveVars isa Tuple || saveVars isa Vector
        Tuple(saveVars)
    else
        error("Invalid saveVars option: $saveVars")
    end
    
    return NamedTuple{selected_vars}(Tuple(copy(getproperty(array, v)) for v in selected_vars))
end

function FMUSaveFunc(u, t, integrator, saveVars)
    input_func!, array, cache, model = integrator.p
    
    input_func!(array, t)
    array.data[1:length(u)] .= u
    EvaluateOutputs!(array, cache, model, t)
    
    return BuildOutputTuple(array, model, saveVars)
end

"""
    simulate(model::AbstractFMUModel, tspan; input_func=(array, t)->nothing, init_func=(array, model)->nothing, solver=Tsit5(), start_from_trim::Bool=false, saveVars=:summary)

Run a native Julia SciML simulation, returning a zero-overhead StructArray of outputs.
"""
function simulate(model::AbstractFMUModel{T}, tspan; 
                  input_func=(array, t)->nothing, 
                  init_func=(array, model)->nothing, 
                  solver=Tsit5(), 
                  start_from_trim::Bool=false, 
                  saveVars=:summary) where T
    array = CreateFMUArray(model)
    cache = AllocateFMUCaches(model, T)
    
    input_func(array, 0.0)
    init_func(array, model)
    InitializeFMU!(array, cache, model, 0.0; start_from_trim=start_from_trim)
    
    n_states = NumberOfStates(model)
    u0 = copy(view(array.data, 1:n_states))
    
    p = (input_func, array, cache, model)
    
    function ode_func!(du, u, p, t)
        input_func, array, cache, model = p
        
        array.data[1:length(u)] .= u
        input_func(array, t)
        EvaluateDerivatives!(du, array, cache, model, t)
    end
    
    prob = ODEProblem(ode_func!, u0, tspan, p)
    
    save_cb_func = (u, t, integrator) -> FMUSaveFunc(u, t, integrator, saveVars)
    dummy_out = save_cb_func(u0, 0.0, (p=p, f=ode_func!))
    outputs = SavedValues(Float64, typeof(dummy_out))
    scb = SavingCallback(save_cb_func, outputs)
    
    sol = solve(prob, solver; callback=scb)
    
    return (sol=sol, outputs=StructArray(outputs.saveval), t=outputs.t)
end
