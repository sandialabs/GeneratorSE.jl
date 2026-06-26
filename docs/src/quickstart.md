# Quick Start

The detailed generator functions expose many design parameters. For integration
work, start by loading the package and inspecting the function docstrings and
examples for the model family you need.

```julia
using GeneratorSE

# Radial-flux PMSG sizing model.
# See the PMSG_arms docstring for required geometric, electrical, and material inputs.
# result = PMSG_arms(...)

# Dynamic model parameters and stepping utilities are exported as well.
# params = PMSG_dynamics_params(...)
# dyn = PMSG_dynamics(...)
```

SIRENOpt wraps this lower-level API in `GeneratorSEModel`,
`generatorse_pmsg_arms_model`, and `generatorse_output_kw` so the system model can
swap generator fidelity without changing controller or storage code.
