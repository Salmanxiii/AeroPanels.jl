# fmi3PlatformTypes.jl
module fmi3PlatformTypes

export fmi3Instance, fmi3InstanceEnvironment, fmi3FMUState, fmi3ValueReference
export fmi3Float32, fmi3Float64, fmi3Int8, fmi3UInt8, fmi3Int16, fmi3UInt16
export fmi3Int32, fmi3UInt32, fmi3Int64, fmi3UInt64, fmi3Boolean, fmi3Char
export fmi3String, fmi3Byte, fmi3Binary, fmi3Clock
export fmi3True, fmi3False, fmi3ClockActive, fmi3ClockInactive

# Core Pointers
const fmi3Instance            = Ptr{Cvoid}
const fmi3InstanceEnvironment = Ptr{Cvoid}
const fmi3FMUState            = Ptr{Cvoid}
const fmi3ValueReference      = Cuint

# Variable Types
const fmi3Float32 = Cfloat
const fmi3Float64 = Cdouble
const fmi3Int8    = Int8
const fmi3UInt8   = UInt8
const fmi3Int16   = Int16
const fmi3UInt16  = UInt16
const fmi3Int32   = Cint
const fmi3UInt32  = Cuint
const fmi3Int64   = Clonglong
const fmi3UInt64  = Culonglong
const fmi3Boolean = Bool
const fmi3Char    = Cchar
const fmi3String  = Ptr{Cchar}
const fmi3Byte    = UInt8
const fmi3Binary  = Ptr{Cvoid}
const fmi3Clock   = Bool

# Constants
const fmi3True          = true
const fmi3False         = false
const fmi3ClockActive   = true
const fmi3ClockInactive = false

end