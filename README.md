# GeneratorSE.jl

[![Tests](https://github.com/sandialabs/GeneratorSE.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/sandialabs/GeneratorSE.jl/actions/workflows/CI.yml)
[![Docs](https://github.com/sandialabs/GeneratorSE.jl/actions/workflows/Docs.yml/badge.svg?branch=master)](https://sandialabs.github.io/GeneratorSE.jl/dev/)
[![Coverage](https://codecov.io/gh/sandialabs/GeneratorSE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sandialabs/GeneratorSE.jl)

Documentation is hosted at [sandialabs.github.io/GeneratorSE.jl/dev/](https://sandialabs.github.io/GeneratorSE.jl/dev/).

GeneratorSE is a set of analytical tools for sizing variable speed wind turbine Generators. The analytical framework involves electromagnetic, structural, and basic thermal design that are integrated to provide the optimal generator design dimensions using conventional magnetic circuit laws.

The tool mainly considers available torque, mechanical power, normal and shear stresses, material properties, and costs to optimize designs of variable-speed wind turbine generators by satisfying specific design criteria.

Original Author: [NREL WISDEM Team](mailto:systems.engineering@nrel.gov)

Julia Changes Author: Kevin Moore, Sandia National Labs

## Documentation

See the [hosted documentation](https://sandialabs.github.io/GeneratorSE.jl/dev/) or the local documentation in the `docs` directory. Note that this translation does not include the interface to OpenMDAO, though this could be restored in the future via [OpenMDAO.jl](https://github.com/byuflowlab/OpenMDAO.jl).

## Modeling Deviations From Original Code

The Julia `PMSG_arms` and `PMSG_axial` models intentionally deviate from the
original translated electrical equations in two places. First, `p` is treated as
pole pairs throughout the models, so electrical angular speed is computed as
`omega_e = p * omega_m` rather than `p * omega_m / 2`. Second, the stator-current
reactive term divides by electrical reactance, `omega_e * L_s`; the original
translated expression effectively divided by `(omega_e * L_s)^2` before
squaring the current component. These changes are dimensional-consistency
corrections and affect calculated current, electrical loading, copper loss, and
efficiency. The Halbach axial model uses the same corrected current expression.

## Installation

`GeneratorSE.jl` is distributed as an unregistered Julia package. Install it
from the public repository URL:

```julia
using Pkg
Pkg.add(url = "https://github.com/sandialabs/GeneratorSE.jl")
```

For local development from a checkout:

```julia
using Pkg
Pkg.develop(path = "/path/to/GeneratorSE.jl")
```

Testing can be done via:

```julia
using Pkg
Pkg.test("GeneratorSE")
```

Legacy optimization scripts live in `analysis/`. They are not part of the
package test target and may require additional solver, plotting, and data-file
dependencies.

## Public API

The main exported sizing and dynamics functions are `PMSG_arms`,
`PMSG_axial`, `PMSG_axial_Halbach`, `PMSG_dynamics`,
`PMSG_dynamics_params`, `segmented_halbach_winding_properties`,
`coreless_winding_inductance`, `copper_resistivity_at_temperature`,
`generator_costing_complex`, and
`generator_costing_simple`.
