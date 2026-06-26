# Theory

GeneratorSE uses analytical electromagnetic, structural, thermal, and costing
relationships to size permanent magnet synchronous generators. The models combine
magnetic circuit assumptions with stress and material constraints to estimate
geometry, losses, mass, and cost from the requested operating point.

For coupled dynamic-platform studies, the package should be treated as a generator
subsystem rather than a complete electrical plant. Rotor torque and speed enter
from the aerodynamic or hydrokinetic model; electrical power, losses, and design
properties pass to converter, battery, and cost models.

The present Julia package does not include the historical OpenMDAO interface.
That is intentional for SIRENOpt integration: the ontology consumes plain Julia
functions so automatic differentiation and package-level testing remain direct.
