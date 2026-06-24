# fmi2PlatformTypes.jl
module fmi2PlatformTypes

export fmi2Component, fmi2ComponentEnvironment, fmi2FMUstate, fmi2ValueReference
export fmi2Real, fmi2Integer, fmi2Boolean, fmi2String
export fmi2True, fmi2False
export fmi2Status, fmi2OK, fmi2Warning, fmi2Discard, fmi2Error, fmi2Fatal, fmi2Pending

# Core Pointers
const fmi2Component            = Ptr{Cvoid}
const fmi2ComponentEnvironment = Ptr{Cvoid}
const fmi2FMUstate             = Ptr{Cvoid}
const fmi2ValueReference       = Cuint

# Variable Types
const fmi2Real                 = Cdouble
const fmi2Integer              = Int32
const fmi2Boolean              = Int32
const fmi2String               = Ptr{Cchar}
const fmi2Status               = Cint

# Constants
const fmi2True                 = Int32(1)
const fmi2False                = Int32(0)

const fmi2OK                   = Cint(0)
const fmi2Warning              = Cint(1)
const fmi2Discard              = Cint(2)
const fmi2Error                = Cint(3)
const fmi2Fatal                = Cint(4)
const fmi2Pending              = Cint(5)

end