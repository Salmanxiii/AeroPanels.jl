module {{FMU_NAME}}

include("fmi3PlatformTypes.jl")
include("fmi3FunctionTypes.jl")

using .fmi3PlatformTypes
using .fmi3FunctionTypes

using LinearAlgebra
using StaticArrays
using AeroPanels

# --- INJECTED CODE START ---
{{USER_CODE}}
# --- INJECTED CODE END ---

"""
A stable instance wrapper for the FMU to ensure perfect type inference for JuliaC.
"""
struct ModelInstance
    model::UnsteadyAeroModel2D{Float64}
    cache::UnsteadyAeroCache{Float64}
    inputs::AeroInputs{Float64}
end

const INSTANCES = Dict{UInt64, ModelInstance}()
const INSTANCE_LOCK = ReentrantLock()
const NEXT_ID = Ref{UInt64}(1)

function get_instance(ptr::fmi3Instance)::ModelInstance
    return INSTANCES[UInt64(ptr)]
end

# Status codes
const fmi3OK = fmi3Status(0)
const fmi3Warning = fmi3Status(1)
const fmi3Discard = fmi3Status(2)
const fmi3Error = fmi3Status(3)
const fmi3Fatal = fmi3Status(4)

# FMI 3.0 Mandatory Functions
Base.@ccallable function fmi3GetVersion()::fmi3String
    return pointer("3.0\0")
end

Base.@ccallable function fmi3SetDebugLogging(instance::fmi3Instance, loggingOn::fmi3Boolean, nCategories::Csize_t, categories::Ptr{fmi3String})::fmi3Status
    return fmi3OK
end

Base.@ccallable function fmi3InstantiateModelExchange(
    instanceName::fmi3String, instantiationToken::fmi3String, resourcePath::fmi3String,
    visible::fmi3Boolean, loggingOn::fmi3Boolean, instanceEnvironment::fmi3InstanceEnvironment,
    logMessage::Ptr{Cvoid}
)::fmi3Instance
    lock(INSTANCE_LOCK) do
        id = NEXT_ID[]
        NEXT_ID[] += 1
        model = {{BUILDER_NAME}}()
        cache = CreateCacheArrays(model)
        inputs = AeroInputs(model)
        INSTANCES[id] = ModelInstance(model, cache, inputs)
        return fmi3Instance(id)
    end
end

Base.@ccallable function fmi3InstantiateCoSimulation(
    instanceName::fmi3String, instantiationToken::fmi3String, resourcePath::fmi3String,
    visible::fmi3Boolean, loggingOn::fmi3Boolean, eventModeUsed::fmi3Boolean,
    earlyReturnAllowed::fmi3Boolean, requiredIntermediateVariables::Ptr{fmi3ValueReference},
    nRequiredIntermediateVariables::Csize_t, instanceEnvironment::fmi3InstanceEnvironment,
    logMessage::Ptr{Cvoid}, intermediateUpdate::Ptr{Cvoid}
)::fmi3Instance
    return fmi3InstantiateModelExchange(instanceName, instantiationToken, resourcePath, visible, loggingOn, instanceEnvironment, logMessage)
end

Base.@ccallable function fmi3FreeInstance(instance::fmi3Instance)::Cvoid
    lock(INSTANCE_LOCK) do
        delete!(INSTANCES, UInt64(instance))
    end
    return nothing
end

Base.@ccallable function fmi3EnterInitializationMode(instance::fmi3Instance, toleranceDefined::fmi3Boolean, tolerance::fmi3Float64, startTime::fmi3Float64, stopTimeDefined::fmi3Boolean, stopTime::fmi3Float64)::fmi3Status; return fmi3OK; end
Base.@ccallable function fmi3ExitInitializationMode(instance::fmi3Instance)::fmi3Status; return fmi3OK; end
Base.@ccallable function fmi3EnterEventMode(instance::fmi3Instance)::fmi3Status; return fmi3OK; end
Base.@ccallable function fmi3Terminate(instance::fmi3Instance)::fmi3Status; return fmi3OK; end
Base.@ccallable function fmi3Reset(instance::fmi3Instance)::fmi3Status
    comp = get_instance(instance)
    fill!(comp.cache.Γw1, 0.0)
    # Reset inputs to default
    comp.inputs.ρ = 1.225
    return fmi3OK
end

Base.@ccallable function fmi3GetFloat64(instance::fmi3Instance, vr::Ptr{fmi3ValueReference}, nvr::Csize_t, values::Ptr{fmi3Float64}, nValues::Csize_t)::fmi3Status
    comp = get_instance(instance)
    vrs = unsafe_wrap(Array, vr, nvr, own=false)
    vals = unsafe_wrap(Array, values, nValues, own=false)

    idx = 1
    for i in 1:nvr
        v = vrs[i]
        if v == 1 # x (States: Gamma_w1)
            res = comp.cache.Γw1
            for j in 1:length(res); vals[idx] = res[j]; idx += 1; end
        elseif v == 2 # u (Inputs)
            # Not typically called for inputs, but provided for completeness
            # Inputs are stored in Body Frame in comp.inputs
            vals[idx] = comp.inputs.vb[1]; idx += 1
            vals[idx] = comp.inputs.vb[2]; idx += 1
            vals[idx] = comp.inputs.vb[3]; idx += 1
            # ... skipping the rest for brevity as it's rare
        elseif v == 3 # y (Outputs: Forces & Moments)
            # Ensure cache is updated with current inputs and states
            # comp.inputs contains Body Frame values
            AeroSolve!(comp.cache, SVector(comp.inputs.vb), SVector(comp.inputs.ab), 
                       SVector(comp.inputs.ωb), SVector(comp.inputs.dωb), 
                       comp.inputs.δc, comp.inputs.dδc, comp.inputs.ddδc, 
                       comp.inputs.rs, comp.inputs.vs, comp.inputs.as, comp.cache.Γw1, comp.model, comp.inputs.ρ)
            
            # GetTotalForces returns results in Geometry Frame
            Fg, Mg = GetTotalForces(comp.cache, comp.model)
            
            # Use AeroPanels helper to convert Geometry to Body Frame
            Fb = AeroPanels.GeometryToBodyAxis(Fg, comp.model.modelProperties)
            Mb = AeroPanels.GeometryToBodyAxis(Mg, comp.model.modelProperties)
            
            vals[idx] = Fb[1]; idx += 1
            vals[idx] = Fb[2]; idx += 1
            vals[idx] = Fb[3]; idx += 1
            vals[idx] = Mb[1]; idx += 1
            vals[idx] = Mb[2]; idx += 1
            vals[idx] = Mb[3]; idx += 1
        elseif v == 4 # der(x)
            # Derivative is in Geometry Frame (consistent with Gamma_w1)
            res = comp.cache.dΓw1
            for j in 1:length(res); vals[idx] = res[j]; idx += 1; end
        else
            return fmi3Error
        end
    end
    return fmi3OK
end

Base.@ccallable function fmi3SetFloat64(instance::fmi3Instance, vr::Ptr{fmi3ValueReference}, nvr::Csize_t, values::Ptr{fmi3Float64}, nValues::Csize_t)::fmi3Status
    comp = get_instance(instance)
    vrs = unsafe_wrap(Array, vr, nvr, own=false)
    vals = unsafe_wrap(Array, values, nValues, own=false)

    idx = 1
    for i in 1:nvr
        v = vrs[i]
        if v == 2 # u: [vx, vy, vz, p, q, r, ax, ay, az, dp, dq, dr, rho]
            # We treat FMI inputs as being in the Body Frame
            comp.inputs.vb .= SVector(vals[idx], vals[idx+1], vals[idx+2]); idx += 3
            comp.inputs.ωb .= SVector(vals[idx], vals[idx+1], vals[idx+2]); idx += 3
            comp.inputs.ab .= SVector(vals[idx], vals[idx+1], vals[idx+2]); idx += 3
            comp.inputs.dωb .= SVector(vals[idx], vals[idx+1], vals[idx+2]); idx += 3
            comp.inputs.ρ = vals[idx]; idx += 1
        else
            return fmi3Error
        end
    end
    return fmi3OK
end

Base.@ccallable function fmi3SetTime(instance::fmi3Instance, time::fmi3Float64)::fmi3Status
    return fmi3OK
end

Base.@ccallable function fmi3SetContinuousStates(instance::fmi3Instance, continuousStates::Ptr{fmi3Float64}, nContinuousStates::Csize_t)::fmi3Status
    comp = get_instance(instance)
    states = unsafe_wrap(Array, continuousStates, nContinuousStates, own=false)
    copyto!(comp.cache.Γw1, states)
    return fmi3OK
end

Base.@ccallable function fmi3GetContinuousStateDerivatives(instance::fmi3Instance, derivatives::Ptr{fmi3Float64}, nContinuousStates::Csize_t)::fmi3Status
    comp = get_instance(instance)
    ders = unsafe_wrap(Array, derivatives, nContinuousStates, own=false)
    
    AeroSolve!(comp.cache, SVector(comp.inputs.vb), SVector(comp.inputs.ab), 
               SVector(comp.inputs.ωb), SVector(comp.inputs.dωb), 
               comp.inputs.δc, comp.inputs.dδc, comp.inputs.ddδc, 
               comp.inputs.rs, comp.inputs.vs, comp.inputs.as, comp.cache.Γw1, comp.model, comp.inputs.ρ)
               
    copyto!(ders, comp.cache.dΓw1)
    return fmi3OK
end

Base.@ccallable function fmi3GetContinuousStates(instance::fmi3Instance, continuousStates::Ptr{fmi3Float64}, nContinuousStates::Csize_t)::fmi3Status
    comp = get_instance(instance)
    states = unsafe_wrap(Array, continuousStates, nContinuousStates, own=false)
    copyto!(states, comp.cache.Γw1)
    return fmi3OK
end

Base.@ccallable function fmi3GetNominalsOfContinuousStates(instance::fmi3Instance, nominals::Ptr{fmi3Float64}, nContinuousStates::Csize_t)::fmi3Status
    noms = unsafe_wrap(Array, nominals, nContinuousStates, own=false)
    fill!(noms, 1.0)
    return fmi3OK
end

Base.@ccallable function fmi3GetNumberOfContinuousStates(instance::fmi3Instance, nContinuousStates::Ptr{Csize_t})::fmi3Status
    unsafe_store!(nContinuousStates, Csize_t(length(get_instance(instance).cache.Γw1)))
    return fmi3OK
end

Base.@ccallable function fmi3EnterContinuousTimeMode(instance::fmi3Instance)::fmi3Status; return fmi3OK; end
Base.@ccallable function fmi3CompletedIntegratorStep(instance::fmi3Instance, noSetFMUStatePriorToCurrentPoint::fmi3Boolean, enterEventMode::Ptr{fmi3Boolean}, terminateSimulation::Ptr{fmi3Boolean})::fmi3Status
    unsafe_store!(enterEventMode, fmi3False); unsafe_store!(terminateSimulation, fmi3False); return fmi3OK
end

Base.@ccallable function fmi3UpdateDiscreteStates(instance::fmi3Instance, discreteStatesNeedUpdate::Ptr{fmi3Boolean}, terminateSimulation::Ptr{fmi3Boolean}, nominalsOfContinuousStatesChanged::Ptr{fmi3Boolean}, valuesOfContinuousStatesChanged::Ptr{fmi3Boolean}, nextEventTimeDefined::Ptr{fmi3Boolean}, nextEventTime::Ptr{fmi3Float64})::fmi3Status
    unsafe_store!(discreteStatesNeedUpdate, fmi3False); unsafe_store!(terminateSimulation, fmi3False)
    unsafe_store!(nominalsOfContinuousStatesChanged, fmi3False); unsafe_store!(valuesOfContinuousStatesChanged, fmi3False)
    unsafe_store!(nextEventTimeDefined, fmi3False); return fmi3OK
end

Base.@ccallable function fmi3GetNumberOfEventIndicators(instance::fmi3Instance, nEventIndicators::Ptr{Csize_t})::fmi3Status
    unsafe_store!(nEventIndicators, Csize_t(0)); return fmi3OK
end

Base.@ccallable function fmi3DoStep(instance::fmi3Instance, currentCommunicationPoint::fmi3Float64, communicationStepSize::fmi3Float64, noSetFMUStatePriorToCurrentPoint::fmi3Boolean, eventHandlingNeeded::Ptr{fmi3Boolean}, terminateSimulation::Ptr{fmi3Boolean}, earlyReturn::Ptr{fmi3Boolean}, lastSuccessfulTime::Ptr{fmi3Float64})::fmi3Status
    comp = get_instance(instance)

    AeroSolve!(comp.cache, SVector(comp.inputs.vb), SVector(comp.inputs.ab), 
               SVector(comp.inputs.ωb), SVector(comp.inputs.dωb), 
               comp.inputs.δc, comp.inputs.dδc, comp.inputs.ddδc, 
               comp.inputs.rs, comp.inputs.vs, comp.inputs.as, comp.cache.Γw1, comp.model, comp.inputs.ρ)
               
    comp.cache.Γw1 .+= comp.cache.dΓw1 .* communicationStepSize

    unsafe_store!(eventHandlingNeeded, fmi3False)
    unsafe_store!(terminateSimulation, fmi3False)
    unsafe_store!(earlyReturn, fmi3False)
    unsafe_store!(lastSuccessfulTime, currentCommunicationPoint + communicationStepSize)
    return fmi3OK
end

Base.@ccallable function fmi3EnterStepMode(instance::fmi3Instance)::fmi3Status; return fmi3OK; end

end
