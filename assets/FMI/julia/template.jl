module {{FMU_NAME}}

include("fmi2PlatformTypes.jl")

using .fmi2PlatformTypes

using LinearAlgebra
using StaticArrays
using GeometryBasics
using AeroPanels
using OrdinaryDiffEqTsit5
using ForwardDiff
using Serialization

# --- INJECTED CODE START ---
{{USER_CODE}}
# --- INJECTED CODE END ---

{{FMU_LAYOUT_DEF}}

const Dual1 = ForwardDiff.Dual{nothing, Float64, 1}

mutable struct ModelInstance
    model::Any
    FMUarray::Any
    DualFMUarray::Any
    cache::Any
    time::Float64
    outputsSynced::Bool
    derivativesSynced::Bool
    integrator::Any
    integrator_dirty::Bool
end

const INSTANCES = Dict{UInt64, ModelInstance}()
const INSTANCE_LOCK = ReentrantLock()
const NEXT_ID = Ref{UInt64}(1)

struct ModelStateSnapshot
    data::Vector{Float64}
    time::Float64
end

const SAVED_STATES = Dict{UInt64, ModelStateSnapshot}()
const STATE_LOCK = ReentrantLock()
const NEXT_STATE_ID = Ref{UInt64}(1)

function get_instance(ptr::fmi2Component)::ModelInstance
    return INSTANCES[UInt64(ptr)]
end

# FMI 2.0 Mandatory Functions
Base.@ccallable function fmi2GetTypesPlatform()::fmi2String
    return pointer("default\0")
end

Base.@ccallable function fmi2GetVersion()::fmi2String
    return pointer("2.0\0")
end

Base.@ccallable function fmi2SetDebugLogging(instance::fmi2Component, loggingOn::fmi2Boolean, nCategories::Csize_t, categories::Ptr{fmi2String})::fmi2Status
    return fmi2OK
end

function ode_wrapper!(du, u, params, t)
    model, FMUarray, cache = params
    nStates = NumberOfStates(model)
    if nStates > 0
        copyto!(view(FMUarray.data, 1:nStates), u)
    end
    EvaluateDerivatives!(du, FMUarray, cache, model, t)
    return nothing
end

Base.@ccallable function fmi2Instantiate(
    instanceName::fmi2String, fmuType::Cint, fmuGUID::fmi2String,
    fmuResourceLocation::fmi2String, callbacks::Ptr{Cvoid},
    visible::fmi2Boolean, loggingOn::fmi2Boolean
)::fmi2Component
    res_dir = ""
    if fmuResourceLocation != Ptr{UInt8}(C_NULL)
        loc_str = unsafe_string(fmuResourceLocation)
        if startswith(loc_str, "file:///")
            res_dir = loc_str[9:end]
        elseif startswith(loc_str, "file://")
            res_dir = loc_str[8:end]
        else
            res_dir = loc_str
        end
        res_dir = replace(res_dir, "/" => "\\")
    end

    model = create_model(res_dir)
    FMUarray = CreateFMUArray(model)
    
    # Pre-allocate Dual FMUarray for directional derivatives
    dual_data = Vector{Dual1}(undef, length(FMUarray.data))
    for i in 1:length(FMUarray.data)
        dual_data[i] = Dual1(FMUarray.data[i], ForwardDiff.Partials((0.0,)))
    end
    DualFMUarray = CreateFMUArray(model, dual_data)
    
    cache = AllocateFMUCaches(model, Float64)
    
    nStates = NumberOfStates(model)
    u0 = nStates > 0 ? Vector{Float64}(view(FMUarray.data, 1:nStates)) : Float64[]
    
    prob = ODEProblem(ode_wrapper!, u0, (0.0, 1e6), (model, FMUarray, cache))
    integrator = init(prob, Tsit5(); reltol=1e-6, abstol=1e-6, save_everystep=false)
    
    inst = ModelInstance(
        model,
        FMUarray,
        DualFMUarray,
        cache,
        0.0,
        false,
        false,
        integrator,
        true
    )
    
    id = lock(INSTANCE_LOCK) do
        cur = NEXT_ID[]
        NEXT_ID[] += 1
        return cur
    end
    
    ptr = Ptr{Cvoid}(id)
    INSTANCES[id] = inst
    return ptr
end

Base.@ccallable function fmi2FreeInstance(instance::fmi2Component)::Cvoid
    if instance != C_NULL
        id = UInt64(instance)
        lock(INSTANCE_LOCK) do
            delete!(INSTANCES, id)
        end
    end
    return nothing
end

Base.@ccallable function fmi2SetupExperiment(
    instance::fmi2Component, toleranceDefined::fmi2Boolean,
    tolerance::fmi2Real, startTime::fmi2Real,
    stopTimeDefined::fmi2Boolean, stopTime::fmi2Real
)::fmi2Status
    comp = get_instance(instance)
    comp.time = startTime
    comp.integrator_dirty = true
    return fmi2OK
end

Base.@ccallable function fmi2EnterInitializationMode(instance::fmi2Component)::fmi2Status
    return fmi2OK
end

Base.@ccallable function fmi2ExitInitializationMode(instance::fmi2Component)::fmi2Status
    comp = get_instance(instance)
    InitializeFMU!(comp.FMUarray, comp.cache, comp.model, comp.time)
    comp.outputsSynced = false
    comp.derivativesSynced = false
    return fmi2OK
end

Base.@ccallable function fmi2Terminate(instance::fmi2Component)::fmi2Status
    return fmi2OK
end

Base.@ccallable function fmi2Reset(instance::fmi2Component)::fmi2Status
    comp = get_instance(instance)
    fill!(comp.FMUarray.data, 0.0)
    comp.time = 0.0
    comp.outputsSynced = false
    comp.derivativesSynced = false
    comp.integrator_dirty = true
    return fmi2OK
end

Base.@ccallable function fmi2GetReal(instance::fmi2Component, vr::Ptr{fmi2ValueReference}, nvr::Csize_t, value::Ptr{fmi2Real})::fmi2Status
    if instance == C_NULL || vr == C_NULL || value == C_NULL
        return fmi2Error
    end
    comp = get_instance(instance)
    vrs = unsafe_wrap(Array, vr, nvr, own=false)
    vals = unsafe_wrap(Array, value, nvr, own=false)
    layout = GetFMULayout(comp.model)
    nStates = NumberOfStates(comp.model)
    out_start = layout.OutputsStartIndex

    any_outputs = false
    any_derivatives = false
    for i in 1:nvr
        v = vrs[i]
        if v >= out_start
            any_outputs = true
        elseif v > nStates && v <= 2 * nStates
            any_derivatives = true
        end
    end

    if any_outputs && !comp.outputsSynced
        EvaluateOutputs!(comp.FMUarray, comp.cache, comp.model, comp.time)
        comp.outputsSynced = true
    end
    if any_derivatives && !comp.derivativesSynced
        EvaluateDerivatives!(comp.FMUarray, comp.cache, comp.model, comp.time)
        comp.derivativesSynced = true
    end

    for i in 1:nvr
        v = vrs[i]
        if v >= 1 && v <= length(comp.FMUarray.data)
            vals[i] = comp.FMUarray.data[v]
        else
            return fmi2Error
        end
    end
    return fmi2OK
end

Base.@ccallable function fmi2SetReal(instance::fmi2Component, vr::Ptr{fmi2ValueReference}, nvr::Csize_t, value::Ptr{fmi2Real})::fmi2Status
    if instance == C_NULL || vr == C_NULL || value == C_NULL
        return fmi2Error
    end
    comp = get_instance(instance)
    vrs = unsafe_wrap(Array, vr, nvr, own=false)
    vals = unsafe_wrap(Array, value, nvr, own=false)
    layout = GetFMULayout(comp.model)
    nStates = NumberOfStates(comp.model)
    out_start = layout.OutputsStartIndex

    comp.outputsSynced = false
    comp.derivativesSynced = false
    
    for i in 1:nvr
        v = vrs[i]
        # Read-only guards for derivatives and outputs
        if (v > nStates && v <= 2 * nStates) || v >= out_start
            continue
        elseif v >= 1 && v <= length(comp.FMUarray.data)
            comp.FMUarray.data[v] = vals[i]
        else
            return fmi2Error
        end
    end
    return fmi2OK
end

Base.@ccallable function fmi2SetTime(instance::fmi2Component, time::fmi2Real)::fmi2Status
    comp = get_instance(instance)
    if comp.time != time
        comp.time = time
        comp.outputsSynced = false
        comp.derivativesSynced = false
    end
    return fmi2OK
end

Base.@ccallable function fmi2SetContinuousStates(instance::fmi2Component, continuousStates::Ptr{fmi2Real}, nContinuousStates::Csize_t)::fmi2Status
    comp = get_instance(instance)
    states = unsafe_wrap(Array, continuousStates, nContinuousStates, own=false)
    nStates = NumberOfStates(comp.model)
    if nContinuousStates == nStates
        for i in 1:nStates
            comp.FMUarray.data[i] = states[i]
        end
        comp.outputsSynced = false
        comp.derivativesSynced = false
        comp.integrator_dirty = true
    end
    return fmi2OK
end

Base.@ccallable function fmi2GetDerivatives(instance::fmi2Component, derivatives::Ptr{fmi2Real}, nContinuousStates::Csize_t)::fmi2Status
    comp = get_instance(instance)
    ders = unsafe_wrap(Array, derivatives, nContinuousStates, own=false)
    nStates = NumberOfStates(comp.model)
    
    if !comp.derivativesSynced
        EvaluateDerivatives!(comp.FMUarray, comp.cache, comp.model, comp.time)
        comp.derivativesSynced = true
    end
               
    for i in 1:nStates
        ders[i] = comp.FMUarray.data[nStates + i]
    end
    return fmi2OK
end

Base.@ccallable function fmi2GetContinuousStates(instance::fmi2Component, continuousStates::Ptr{fmi2Real}, nContinuousStates::Csize_t)::fmi2Status
    comp = get_instance(instance)
    states = unsafe_wrap(Array, continuousStates, nContinuousStates, own=false)
    nStates = NumberOfStates(comp.model)
    for i in 1:nStates
        states[i] = comp.FMUarray.data[i]
    end
    return fmi2OK
end

Base.@ccallable function fmi2GetNominalsOfContinuousStates(instance::fmi2Component, nominals::Ptr{fmi2Real}, nContinuousStates::Csize_t)::fmi2Status
    noms = unsafe_wrap(Array, nominals, nContinuousStates, own=false)
    fill!(noms, 1.0)
    return fmi2OK
end

Base.@ccallable function fmi2EnterContinuousTimeMode(instance::fmi2Component)::fmi2Status; return fmi2OK; end
Base.@ccallable function fmi2CompletedIntegratorStep(instance::fmi2Component, noSetFMUStatePriorToCurrentPoint::fmi2Boolean, enterEventMode::Ptr{fmi2Boolean}, terminateSimulation::Ptr{fmi2Boolean})::fmi2Status
    unsafe_store!(enterEventMode, fmi2False); unsafe_store!(terminateSimulation, fmi2False); return fmi2OK
end

Base.@ccallable function fmi2NewDiscreteStates(instance::fmi2Component, eventInfo::Ptr{Cvoid})::fmi2Status
    return fmi2OK
end

Base.@ccallable function fmi2GetFMUstate(instance::fmi2Component, FMUstate::Ptr{Ptr{Cvoid}})::fmi2Status
    if instance == C_NULL || FMUstate == C_NULL
        return fmi2Error
    end
    comp = get_instance(instance)
    snapshot = ModelStateSnapshot(copy(comp.FMUarray.data), comp.time)
    
    id = lock(STATE_LOCK) do
        cur = NEXT_STATE_ID[]
        NEXT_STATE_ID[] += 1
        SAVED_STATES[cur] = snapshot
        return cur
    end
    
    unsafe_store!(FMUstate, Ptr{Cvoid}(id))
    return fmi2OK
end

Base.@ccallable function fmi2SetFMUstate(instance::fmi2Component, FMUstate::Ptr{Cvoid})::fmi2Status
    if instance == C_NULL || FMUstate == C_NULL
        return fmi2Error
    end
    comp = get_instance(instance)
    id = UInt64(FMUstate)
    snapshot = lock(STATE_LOCK) do
        get(SAVED_STATES, id, nothing)
    end
    if snapshot === nothing
        return fmi2Error
    end
    copyto!(comp.FMUarray.data, snapshot.data)
    comp.time = snapshot.time
    comp.outputsSynced = false
    comp.derivativesSynced = false
    comp.integrator_dirty = true
    return fmi2OK
end

Base.@ccallable function fmi2FreeFMUstate(instance::fmi2Component, FMUstate::Ptr{Ptr{Cvoid}})::fmi2Status
    if FMUstate == C_NULL
        return fmi2Error
    end
    state_ptr = unsafe_load(FMUstate)
    if state_ptr == C_NULL
        return fmi2OK
    end
    id = UInt64(state_ptr)
    lock(STATE_LOCK) do
        delete!(SAVED_STATES, id)
    end
    unsafe_store!(FMUstate, Ptr{Cvoid}(0))
    return fmi2OK
end

# Stubs for unused types
Base.@ccallable function fmi2GetInteger(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{Int32})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetBoolean(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2Boolean})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetString(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2String})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2SetInteger(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{Int32})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2SetBoolean(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2Boolean})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2SetString(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2String})::fmi2Status; return fmi2Error; end

Base.@ccallable function fmi2SerializedFMUstateSize(p1::fmi2Component, p2::Ptr{Cvoid}, p3::Ptr{Csize_t})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2SerializeFMUstate(p1::fmi2Component, p2::Ptr{Cvoid}, p3::Ptr{UInt8}, p4::Csize_t)::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2DeSerializeFMUstate(p1::fmi2Component, p2::Ptr{UInt8}, p3::Csize_t, p4::Ptr{Ptr{Cvoid}})::fmi2Status; return fmi2Error; end

Base.@ccallable function fmi2GetDirectionalDerivative(
    instance::fmi2Component,
    vUnknown::Ptr{fmi2ValueReference}, nUnknown::Csize_t,
    vKnown::Ptr{fmi2ValueReference}, nKnown::Csize_t,
    dvKnown::Ptr{fmi2Real},
    dvUnknown::Ptr{fmi2Real}
)::fmi2Status
    if instance == C_NULL || vUnknown == C_NULL || vKnown == C_NULL || dvKnown == C_NULL || dvUnknown == C_NULL
        return fmi2Error
    end
    comp = get_instance(instance)
    
    vrs_known = unsafe_wrap(Array, vKnown, nKnown, own=false)
    seeds = unsafe_wrap(Array, dvKnown, nKnown, own=false)
    vrs_unknown = unsafe_wrap(Array, vUnknown, nUnknown, own=false)
    outputs = unsafe_wrap(Array, dvUnknown, nUnknown, own=false)
    
    if all(s -> s == 0.0, seeds)
        fill!(outputs, 0.0)
        return fmi2OK
    end

    # Copy nominal Float64 values into DualFMUarray before seeding
    for i in 1:length(comp.FMUarray.data)
        comp.DualFMUarray.data[i] = Dual1(comp.FMUarray.data[i], ForwardDiff.Partials((0.0,)))
    end

    GetDirectionalDerivative!(
        outputs,
        vrs_unknown,
        vrs_known,
        seeds,
        comp.DualFMUarray,
        comp.model,
        comp.cache,
        comp.time
    )

    return fmi2OK
end

# Co-Simulation Functions
Base.@ccallable function fmi2SetRealInputDerivatives(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{Int32}, p5::Ptr{fmi2Real})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetRealOutputDerivatives(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{Int32}, p5::Ptr{fmi2Real})::fmi2Status; return fmi2Error; end

Base.@ccallable function fmi2DoStep(
    instance::fmi2Component,
    currentCommunicationPoint::fmi2Real,
    communicationStepSize::fmi2Real,
    noSetFMUStatePriorToCurrentPoint::fmi2Boolean
)::fmi2Status
    comp = get_instance(instance)
    nStates = NumberOfStates(comp.model)

    if comp.integrator_dirty || comp.integrator.t != currentCommunicationPoint
        if nStates > 0
            u0 = Vector{Float64}(view(comp.FMUarray.data, 1:nStates))
            reinit!(comp.integrator, u0; t0=currentCommunicationPoint, erase_sol=true, reset_dt=true)
        end
        comp.integrator_dirty = false
    end
    
    if nStates > 0
        step!(comp.integrator, communicationStepSize, true)
        copyto!(view(comp.FMUarray.data, 1:nStates), comp.integrator.u)
    end
    
    comp.time = currentCommunicationPoint + communicationStepSize
    comp.outputsSynced = false
    comp.derivativesSynced = false
    return fmi2OK
end

Base.@ccallable function fmi2CancelStep(p1::fmi2Component)::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetStatus(p1::fmi2Component, p2::Cint, p3::Ptr{Cint})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetRealStatus(p1::fmi2Component, p2::Cint, p3::Ptr{fmi2Real})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetIntegerStatus(p1::fmi2Component, p2::Cint, p3::Ptr{Int32})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetBooleanStatus(p1::fmi2Component, p2::Cint, p3::Ptr{fmi2Boolean})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetStringStatus(p1::fmi2Component, p2::Cint, p3::Ptr{fmi2String})::fmi2Status; return fmi2Error; end

end
