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

## Electrical Equation Conventions

The axial and arms PMSG functions use `p` as pole pairs. Electrical angular
speed is therefore

```math
\omega_e = p \omega_m
```

rather than `p * omega_m / 2`. This is an intentional deviation from the
original translated code, which undercounted electrical angular speed when `p`
already represented pole pairs.

The steady-state current estimate also treats `omega_e * L_s` as electrical
reactance. The corrected reactive-current component is proportional to

```math
\frac{E_p - \sqrt{G}}{\omega_e L_s}
```

before forming the current magnitude. The earlier translated expression divided
by `(omega_e * L_s)^2` inside that component. The corrected form is dimensionally
consistent and changes predicted current, electrical loading, losses, and
efficiency for `PMSG_arms`, `PMSG_axial`, and the shared Halbach axial current
calculation.
