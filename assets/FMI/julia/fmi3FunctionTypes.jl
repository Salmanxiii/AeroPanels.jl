# fmi3FunctionTypes.jl
module fmi3FunctionTypes

export fmi3Status, fmi3OK, fmi3Warning, fmi3Discard, fmi3Error, fmi3Fatal
export fmi3DependencyKind, fmi3IntervalQualifier

# fmi3Status Enum mappings
const fmi3Status  = Cint
const fmi3OK      = Cint(0)
const fmi3Warning = Cint(1)
const fmi3Discard = Cint(2)
const fmi3Error   = Cint(3)
const fmi3Fatal   = Cint(4)

# fmi3DependencyKind Enum mappings
const fmi3DependencyKind = Cint
const fmi3Independent    = Cint(0)
const fmi3Constant       = Cint(1)
const fmi3Fixed          = Cint(2)
const fmi3Tunable        = Cint(3)
const fmi3Discrete       = Cint(4)
const fmi3Dependent      = Cint(5)

# fmi3IntervalQualifier Enum mappings
const fmi3IntervalQualifier   = Cint
const fmi3IntervalNotYetKnown = Cint(0)
const fmi3IntervalUnchanged   = Cint(1)
const fmi3IntervalChanged     = Cint(2)

end