# fmi2Functions.jl
module fmi2Functions

include("fmi2PlatformTypes.jl")

using .fmi2PlatformTypes

Base.@ccallable function fmi2GetVersion()::fmi2String
    return pointer("2.0\0")
end

Base.@ccallable function fmi2Instantiate(
    instanceName::fmi2String,
    fmuType::Cint,
    fmuGUID::fmi2String,
    fmuResourceLocation::fmi2String,
    callbacks::Ptr{Cvoid},
    visible::fmi2Boolean,
    loggingOn::fmi2Boolean
)::fmi2Component
    return Ptr{Cvoid}(1) 
end

Base.@ccallable function fmi2FreeInstance(instance::fmi2Component)::Cvoid
    return nothing
end

Base.@ccallable function fmi2EnterInitializationMode(instance::fmi2Component)::fmi2Status
    return fmi2OK
end

Base.@ccallable function fmi2ExitInitializationMode(instance::fmi2Component)::fmi2Status
    return fmi2OK
end

Base.@ccallable function fmi2SetReal(
    instance::fmi2Component,
    valueReferences::Ptr{fmi2ValueReference},
    nValueReferences::Csize_t,
    values::Ptr{fmi2Real}
)::fmi2Status
    return fmi2OK
end

Base.@ccallable function fmi2GetReal(
    instance::fmi2Component,
    valueReferences::Ptr{fmi2ValueReference},
    nValueReferences::Csize_t,
    values::Ptr{fmi2Real}
)::fmi2Status
    return fmi2OK
end

end # module