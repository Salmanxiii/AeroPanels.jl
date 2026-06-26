module {{FMU_NAME}}

include("fmi2PlatformTypes.jl")

using .fmi2PlatformTypes

using LinearAlgebra
using StaticArrays
using GeometryBasics
using AeroPanels
using OrdinaryDiffEqTsit5

# --- INJECTED CODE START ---
{{USER_CODE}}
# --- INJECTED CODE END ---

{{FMU_LAYOUT_DEF}}

"""
A stable instance wrapper for the FMU to ensure perfect type inference for JuliaC.
"""
struct AeroODEParams
    model::UnsteadyAeroModel2D{Float64}
    cache::UnsteadyAeroCache{Float64}
    inputs::AeroInputs{Float64}
end

function aero_ode!(dy, y, params::AeroODEParams, t)
    params.model(dy, y, (params.inputs, params.cache), t)
    copyto!(params.cache.Γw1, y)
    copyto!(params.cache.dΓw1, dy)
    return nothing
end

mutable struct ModelInstance
    model::UnsteadyAeroModel2D{Float64}
    cache::UnsteadyAeroCache{Float64}
    inputs::AeroInputs{Float64}
    time::Float64
    outputsSynced::Bool
    solverSynced::Bool
    cached_Fb::Vector{Float64}
    cached_Mb::Vector{Float64}
    cached_Fmp::Vector{Float64}
    integrator::Any
    integrator_dirty::Bool
end

const INSTANCES = Dict{UInt64, ModelInstance}()
const INSTANCE_LOCK = ReentrantLock()
const NEXT_ID = Ref{UInt64}(1)

struct ModelStateSnapshot
    Γw1::Vector{Float64}
    vb::Vector{Float64}
    wb::Vector{Float64}
    ab::Vector{Float64}
    dwb::Vector{Float64}
    ρ::Float64
    δc::Vector{Float64}
    dδc::Vector{Float64}
    ddδc::Vector{Float64}
    rs::Vector{Point3{Float64}}
    vs::Vector{Point3{Float64}}
    as::Vector{Point3{Float64}}
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

Base.@ccallable function fmi2Instantiate(
    instanceName::fmi2String, fmuType::Cint, fmuGUID::fmi2String,
    fmuResourceLocation::fmi2String, callbacks::Ptr{Cvoid},
    visible::fmi2Boolean, loggingOn::fmi2Boolean
)::fmi2Component
    lock(INSTANCE_LOCK) do
        id = NEXT_ID[]
        NEXT_ID[] += 1
        model = {{BUILDER_NAME}}()
        cache = CreateCacheArrays(model)
        inputs = AeroInputs(model)
        cached_Fb = zeros(Float64, 3)
        cached_Mb = zeros(Float64, 3)
        cached_Fmp = zeros(Float64, 6 * FMU_LAYOUT.nmp)
        
        ode_params = AeroODEParams(model, cache, inputs)
        prob = ODEProblem(aero_ode!, copy(cache.Γw1), (0.0, 1e9), ode_params)
        integrator = init(prob, Tsit5(); save_everystep=false, adaptive=true, reltol=1e-3, abstol=1e-6)
        
        comp = ModelInstance(model, cache, inputs, 0.0, false, false, cached_Fb, cached_Mb, cached_Fmp, integrator, false)
        
        INSTANCES[id] = comp
        return fmi2Component(id)
    end
end

Base.@ccallable function fmi2FreeInstance(instance::fmi2Component)::Cvoid
    lock(INSTANCE_LOCK) do
        delete!(INSTANCES, UInt64(instance))
    end
    return nothing
end

Base.@ccallable function fmi2SetupExperiment(
    instance::fmi2Component, toleranceDefined::fmi2Boolean, tolerance::fmi2Real,
    startTime::fmi2Real, stopTimeDefined::fmi2Boolean, stopTime::fmi2Real
)::fmi2Status
    return fmi2OK
end

Base.@ccallable function fmi2EnterInitializationMode(instance::fmi2Component)::fmi2Status; return fmi2OK; end
Base.@ccallable function fmi2ExitInitializationMode(instance::fmi2Component)::fmi2Status; return fmi2OK; end
Base.@ccallable function fmi2EnterEventMode(instance::fmi2Component)::fmi2Status; return fmi2OK; end
Base.@ccallable function fmi2Terminate(instance::fmi2Component)::fmi2Status; return fmi2OK; end

Base.@ccallable function fmi2Reset(instance::fmi2Component)::fmi2Status
    comp = get_instance(instance)
    fill!(comp.cache.Γw1, 0.0)
    comp.inputs.ρ = 1.225
    comp.time = 0.0
    comp.outputsSynced = false
    comp.solverSynced = false
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

    any_outputs_requested = false
    for i in 1:nvr
        if vrs[i] >= FMU_LAYOUT.vr_Fb
            any_outputs_requested = true
            break
        end
    end

    if any_outputs_requested && !comp.outputsSynced
        AeroSolve!(comp.cache, SVector(comp.inputs.vb), SVector(comp.inputs.ab), 
                   SVector(comp.inputs.ωb), SVector(comp.inputs.dωb), 
                   comp.inputs.δc, comp.inputs.dδc, comp.inputs.ddδc, 
                   comp.inputs.rs, comp.inputs.vs, comp.inputs.as, comp.cache.Γw1, comp.model, comp.inputs.ρ)
        comp.solverSynced = true
                   
        Fb, Mb, Fmp = AeroPanels.Outputs(comp.cache, comp.model)
        copyto!(comp.cached_Fb, Fb)
        copyto!(comp.cached_Mb, Mb)
        if FMU_LAYOUT.nmp > 0
            copyto!(comp.cached_Fmp, Fmp)
        end
        comp.outputsSynced = true
    end

    for i in 1:nvr
        v = vrs[i]
        
        if v >= 1 && v <= FMU_LAYOUT.nx
            vals[i] = comp.cache.Γw1[v]
            
        elseif v >= FMU_LAYOUT.vr_der_x && v < FMU_LAYOUT.vr_der_x + FMU_LAYOUT.nx
            vals[i] = comp.cache.dΓw1[v - FMU_LAYOUT.vr_der_x + 1]
            
        elseif v >= FMU_LAYOUT.vr_vb && v < FMU_LAYOUT.vr_vb + 3
            vals[i] = comp.inputs.vb[v - FMU_LAYOUT.vr_vb + 1]
            
        elseif v >= FMU_LAYOUT.vr_wb && v < FMU_LAYOUT.vr_wb + 3
            vals[i] = comp.inputs.ωb[v - FMU_LAYOUT.vr_wb + 1]
            
        elseif v >= FMU_LAYOUT.vr_ab && v < FMU_LAYOUT.vr_ab + 3
            vals[i] = comp.inputs.ab[v - FMU_LAYOUT.vr_ab + 1]
            
        elseif v >= FMU_LAYOUT.vr_dwb && v < FMU_LAYOUT.vr_dwb + 3
            vals[i] = comp.inputs.dωb[v - FMU_LAYOUT.vr_dwb + 1]
            
        elseif v == FMU_LAYOUT.vr_rho
            vals[i] = comp.inputs.ρ
            
        elseif v >= FMU_LAYOUT.vr_cs_start && v < FMU_LAYOUT.vr_cs_start + 3 * FMU_LAYOUT.nCtrl
            offset = v - FMU_LAYOUT.vr_cs_start
            cs_idx = div(offset, 3) + 1
            type_idx = mod(offset, 3)
            
            if type_idx == 0
                vals[i] = comp.inputs.δc[cs_idx]
            elseif type_idx == 1
                vals[i] = comp.inputs.dδc[cs_idx]
            else
                vals[i] = comp.inputs.ddδc[cs_idx]
            end
            
        elseif FMU_LAYOUT.nVert > 0 && v >= FMU_LAYOUT.vr_rs && v < FMU_LAYOUT.vr_rs + 3 * FMU_LAYOUT.nVert
            vals[i] = comp.inputs.rs[div(v - FMU_LAYOUT.vr_rs, 3) + 1][mod(v - FMU_LAYOUT.vr_rs, 3) + 1]
            
        elseif FMU_LAYOUT.nVert > 0 && v >= FMU_LAYOUT.vr_vs && v < FMU_LAYOUT.vr_vs + 3 * FMU_LAYOUT.nVert
            vals[i] = comp.inputs.vs[div(v - FMU_LAYOUT.vr_vs, 3) + 1][mod(v - FMU_LAYOUT.vr_vs, 3) + 1]
            
        elseif FMU_LAYOUT.nVert > 0 && v >= FMU_LAYOUT.vr_as && v < FMU_LAYOUT.vr_as + 3 * FMU_LAYOUT.nVert
            vals[i] = comp.inputs.as[div(v - FMU_LAYOUT.vr_as, 3) + 1][mod(v - FMU_LAYOUT.vr_as, 3) + 1]
            
        elseif v >= FMU_LAYOUT.vr_Fb && v < FMU_LAYOUT.vr_Fb + 3
            vals[i] = comp.cached_Fb[v - FMU_LAYOUT.vr_Fb + 1]
            
        elseif v >= FMU_LAYOUT.vr_Mb && v < FMU_LAYOUT.vr_Mb + 3
            vals[i] = comp.cached_Mb[v - FMU_LAYOUT.vr_Mb + 1]
            
        elseif FMU_LAYOUT.nmp > 0 && v >= FMU_LAYOUT.vr_mp_start && v < FMU_LAYOUT.vr_mp_start + 6 * FMU_LAYOUT.nmp
            vals[i] = comp.cached_Fmp[v - FMU_LAYOUT.vr_mp_start + 1]
            
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

    comp.outputsSynced = false
    comp.solverSynced = false
    for i in 1:nvr
        v = vrs[i]
        
        if v >= 1 && v <= FMU_LAYOUT.nx
            comp.cache.Γw1[v] = vals[i]
            
        elseif v >= FMU_LAYOUT.vr_vb && v < FMU_LAYOUT.vr_vb + 3
            comp.inputs.vb[v - FMU_LAYOUT.vr_vb + 1] = vals[i]
            
        elseif v >= FMU_LAYOUT.vr_wb && v < FMU_LAYOUT.vr_wb + 3
            comp.inputs.ωb[v - FMU_LAYOUT.vr_wb + 1] = vals[i]
            
        elseif v >= FMU_LAYOUT.vr_ab && v < FMU_LAYOUT.vr_ab + 3
            comp.inputs.ab[v - FMU_LAYOUT.vr_ab + 1] = vals[i]
            
        elseif v >= FMU_LAYOUT.vr_dwb && v < FMU_LAYOUT.vr_dwb + 3
            comp.inputs.dωb[v - FMU_LAYOUT.vr_dwb + 1] = vals[i]
            
        elseif v == FMU_LAYOUT.vr_rho
            comp.inputs.ρ = vals[i]
            
        elseif v >= FMU_LAYOUT.vr_cs_start && v < FMU_LAYOUT.vr_cs_start + 3 * FMU_LAYOUT.nCtrl
            offset = v - FMU_LAYOUT.vr_cs_start
            cs_idx = div(offset, 3) + 1
            type_idx = mod(offset, 3)
            
            if type_idx == 0
                comp.inputs.δc[cs_idx] = vals[i]
            elseif type_idx == 1
                comp.inputs.dδc[cs_idx] = vals[i]
            else
                comp.inputs.ddδc[cs_idx] = vals[i]
            end
            
        elseif FMU_LAYOUT.nVert > 0 && v >= FMU_LAYOUT.vr_rs && v < FMU_LAYOUT.vr_rs + 3 * FMU_LAYOUT.nVert
            idx = v - FMU_LAYOUT.vr_rs + 1
            node_idx = div(idx - 1, 3) + 1
            coord_idx = mod(idx - 1, 3) + 1
            p = comp.inputs.rs[node_idx]
            coords = [p[1], p[2], p[3]]
            coords[coord_idx] = vals[i]
            comp.inputs.rs[node_idx] = Point3(coords[1], coords[2], coords[3])
            
        elseif FMU_LAYOUT.nVert > 0 && v >= FMU_LAYOUT.vr_vs && v < FMU_LAYOUT.vr_vs + 3 * FMU_LAYOUT.nVert
            idx = v - FMU_LAYOUT.vr_vs + 1
            node_idx = div(idx - 1, 3) + 1
            coord_idx = mod(idx - 1, 3) + 1
            p = comp.inputs.vs[node_idx]
            coords = [p[1], p[2], p[3]]
            coords[coord_idx] = vals[i]
            comp.inputs.vs[node_idx] = Point3(coords[1], coords[2], coords[3])
            
        elseif FMU_LAYOUT.nVert > 0 && v >= FMU_LAYOUT.vr_as && v < FMU_LAYOUT.vr_as + 3 * FMU_LAYOUT.nVert
            idx = v - FMU_LAYOUT.vr_as + 1
            node_idx = div(idx - 1, 3) + 1
            coord_idx = mod(idx - 1, 3) + 1
            p = comp.inputs.as[node_idx]
            coords = [p[1], p[2], p[3]]
            coords[coord_idx] = vals[i]
            comp.inputs.as[node_idx] = Point3(coords[1], coords[2], coords[3])
            
        elseif v >= FMU_LAYOUT.vr_der_x && v < FMU_LAYOUT.vr_der_x + FMU_LAYOUT.nx
            # Silently ignore writes to derivative and outputs to protect against initialization writes
            continue
        elseif v >= FMU_LAYOUT.vr_Fb && v < FMU_LAYOUT.vr_Fb + 3
            continue
        elseif v >= FMU_LAYOUT.vr_Mb && v < FMU_LAYOUT.vr_Mb + 3
            continue
        elseif FMU_LAYOUT.nmp > 0 && v >= FMU_LAYOUT.vr_mp_start && v < FMU_LAYOUT.vr_mp_start + 6 * FMU_LAYOUT.nmp
            continue
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
        comp.solverSynced = false
    end
    return fmi2OK
end

Base.@ccallable function fmi2SetContinuousStates(instance::fmi2Component, continuousStates::Ptr{fmi2Real}, nContinuousStates::Csize_t)::fmi2Status
    comp = get_instance(instance)
    states = unsafe_wrap(Array, continuousStates, nContinuousStates, own=false)
    copyto!(comp.cache.Γw1, states)
    comp.outputsSynced = false
    comp.solverSynced = false
    comp.integrator_dirty = true
    return fmi2OK
end

Base.@ccallable function fmi2GetDerivatives(instance::fmi2Component, derivatives::Ptr{fmi2Real}, nContinuousStates::Csize_t)::fmi2Status
    comp = get_instance(instance)
    ders = unsafe_wrap(Array, derivatives, nContinuousStates, own=false)
    
    if !comp.solverSynced
        comp.model(comp.cache.dΓw1, comp.cache.Γw1, (comp.inputs, comp.cache), comp.time)
        comp.solverSynced = true
    end
               
    copyto!(ders, comp.cache.dΓw1)
    return fmi2OK
end

Base.@ccallable function fmi2GetContinuousStates(instance::fmi2Component, continuousStates::Ptr{fmi2Real}, nContinuousStates::Csize_t)::fmi2Status
    comp = get_instance(instance)
    states = unsafe_wrap(Array, continuousStates, nContinuousStates, own=false)
    copyto!(states, comp.cache.Γw1)
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
    if eventInfo == C_NULL
        return fmi2Error
    end
    unsafe_store!(convert(Ptr{Int32}, eventInfo + 0), Int32(0))  # newDiscreteStatesNeeded = false
    unsafe_store!(convert(Ptr{Int32}, eventInfo + 4), Int32(0))  # terminateSimulation = false
    unsafe_store!(convert(Ptr{Int32}, eventInfo + 8), Int32(0))  # nominalsOfContinuousStatesChanged = false
    unsafe_store!(convert(Ptr{Int32}, eventInfo + 12), Int32(0)) # valuesOfContinuousStatesChanged = false
    unsafe_store!(convert(Ptr{Int32}, eventInfo + 16), Int32(0)) # nextEventTimeDefined = false
    return fmi2OK
end

Base.@ccallable function fmi2GetEventIndicators(instance::fmi2Component, eventIndicators::Ptr{fmi2Real}, ni::Csize_t)::fmi2Status
    return fmi2OK
end

# Stubs for state getting/setting and serialization
Base.@ccallable function fmi2GetFMUstate(instance::fmi2Component, FMUstate::Ptr{Ptr{Cvoid}})::fmi2Status
    if instance == C_NULL || FMUstate == C_NULL
        return fmi2Error
    end
    comp = get_instance(instance)
    snapshot = ModelStateSnapshot(
        copy(comp.cache.Γw1),
        copy(comp.inputs.vb),
        copy(comp.inputs.ωb),
        copy(comp.inputs.ab),
        copy(comp.inputs.dωb),
        comp.inputs.ρ,
        copy(comp.inputs.δc),
        copy(comp.inputs.dδc),
        copy(comp.inputs.ddδc),
        copy(comp.inputs.rs),
        copy(comp.inputs.vs),
        copy(comp.inputs.as),
        comp.time
    )
    lock(STATE_LOCK) do
        id = NEXT_STATE_ID[]
        NEXT_STATE_ID[] += 1
        SAVED_STATES[id] = snapshot
        unsafe_store!(FMUstate, Ptr{Cvoid}(id))
    end
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
    copyto!(comp.cache.Γw1, snapshot.Γw1)
    copyto!(comp.inputs.vb, snapshot.vb)
    copyto!(comp.inputs.ωb, snapshot.wb)
    copyto!(comp.inputs.ab, snapshot.ab)
    copyto!(comp.inputs.dωb, snapshot.dwb)
    comp.inputs.ρ = snapshot.ρ
    copyto!(comp.inputs.δc, snapshot.δc)
    copyto!(comp.inputs.dδc, snapshot.dδc)
    copyto!(comp.inputs.ddδc, snapshot.ddδc)
    copyto!(comp.inputs.rs, snapshot.rs)
    copyto!(comp.inputs.vs, snapshot.vs)
    copyto!(comp.inputs.as, snapshot.as)
    comp.time = snapshot.time
    comp.outputsSynced = false
    comp.solverSynced = false
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

# Stubs for other variables and CS
Base.@ccallable function fmi2GetInteger(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{Int32})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetBoolean(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2Boolean})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetString(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2String})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2SetInteger(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{Int32})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2SetBoolean(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2Boolean})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2SetString(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2String})::fmi2Status; return fmi2Error; end

Base.@ccallable function fmi2SerializedFMUstateSize(p1::fmi2Component, p2::Ptr{Cvoid}, p3::Ptr{Csize_t})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2SerializeFMUstate(p1::fmi2Component, p2::Ptr{Cvoid}, p3::Ptr{UInt8}, p4::Csize_t)::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2DeSerializeFMUstate(p1::fmi2Component, p2::Ptr{UInt8}, p3::Csize_t, p4::Ptr{Ptr{Cvoid}})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetDirectionalDerivative(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{fmi2ValueReference}, p5::Csize_t, p6::Ptr{fmi2Real}, p7::Ptr{fmi2Real})::fmi2Status; return fmi2Error; end

# Co-Simulation Stubs
Base.@ccallable function fmi2SetRealInputDerivatives(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{Int32}, p5::Ptr{fmi2Real})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetRealOutputDerivatives(p1::fmi2Component, p2::Ptr{fmi2ValueReference}, p3::Csize_t, p4::Ptr{Int32}, p5::Ptr{fmi2Real})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2DoStep(
    instance::fmi2Component,
    currentCommunicationPoint::fmi2Real,
    communicationStepSize::fmi2Real,
    noSetFMUStatePriorToCurrentPoint::fmi2Boolean
)::fmi2Status
    comp = get_instance(instance)

    if comp.integrator_dirty || comp.integrator.t != currentCommunicationPoint
        reinit!(comp.integrator, comp.cache.Γw1; t0=currentCommunicationPoint, erase_sol=true, reset_dt=true)
        comp.integrator_dirty = false
    end
    
    step!(comp.integrator, communicationStepSize, true)
    
    copyto!(comp.cache.Γw1, comp.integrator.u)
    comp.time = currentCommunicationPoint + communicationStepSize
    comp.outputsSynced = false
    comp.solverSynced = false
    return fmi2OK
end
Base.@ccallable function fmi2CancelStep(p1::fmi2Component)::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetStatus(p1::fmi2Component, p2::Cint, p3::Ptr{Cint})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetRealStatus(p1::fmi2Component, p2::Cint, p3::Ptr{fmi2Real})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetIntegerStatus(p1::fmi2Component, p2::Cint, p3::Ptr{Int32})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetBooleanStatus(p1::fmi2Component, p2::Cint, p3::Ptr{fmi2Boolean})::fmi2Status; return fmi2Error; end
Base.@ccallable function fmi2GetStringStatus(p1::fmi2Component, p2::Cint, p3::Ptr{fmi2String})::fmi2Status; return fmi2Error; end

end
