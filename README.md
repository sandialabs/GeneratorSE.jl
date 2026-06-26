# GeneratorSE.jl

GeneratorSE is a set of analytical tools for sizing variable speed wind turbine Generators. The analytical framework involves electromagnetic, structural, and basic thermal design that are integrated to provide the optimal generator design dimensions using conventional magnetic circuit laws.

The tool mainly considers available torque, mechanical power, normal and shear stresses, material properties, and costs to optimize designs of variable-speed wind turbine generators by satisfying specific design criteria.

Original Author: [NREL WISDEM Team](mailto:systems.engineering@nrel.gov)

Julia Changes Author: Kevin Moore, Sandia National Labs

## Documentation

See local documentation in the `docs`-directory.  Note that this translation does not include the interface to OpenMDAO, though this could be restored in the future via [OpenMDAO.jl](https://github.com/byuflowlab/OpenMDAO.jl)

## Installation

`GeneratorSE.jl` is distributed as an unregistered Julia package. Install it
from the public repository URL:

```julia
using Pkg
Pkg.add(url = "https://github.com/kevmoor/GeneratorSE.jl")
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
`PMSG_dynamics_params`, `generator_costing_complex`, and
`generator_costing_simple`.
