module GeneratorSE
import FLOWMath
path = splitdir(@__FILE__)[1]
export PMSG_arms, PMSG_axial, PMSG_axial_Halbach, PMSG_dynamics, PMSG_dynamics_params, generator_costing_complex, generator_costing_simple#, PMSG_outer

const _GENERATORSE_SMOOTH_DELTA = 1.0e-9

_smooth_abs(x; delta = _GENERATORSE_SMOOTH_DELTA) = FLOWMath.abs_smooth(x, delta)
_smooth_max(a, b; delta = _GENERATORSE_SMOOTH_DELTA) = (a + b + _smooth_abs(a - b; delta)) / 2
_smooth_min(a, b; delta = _GENERATORSE_SMOOTH_DELTA) = (a + b - _smooth_abs(a - b; delta)) / 2

include("$path/PMSG_arms.jl")
include("$path/PMSG_axial.jl")
include("$path/PMSG_axial_Halbach.jl")
include("$path/PMSG_outer.jl")
include("$path/PMSG_dynamics.jl")

end
