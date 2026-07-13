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

## Physical winding and segmented-Halbach options

The generic axial sizing equations retain their historical slot-fill,
distributed-inductance, and harmonic-Halbach defaults. Those defaults are not
appropriate for every hand-wound coreless prototype. `PMSG_axial_Halbach`
therefore has opt-in physical paths that require explicit manufacturing inputs
and do not use measured electrical constants.

With wire area and a turn path, phase resistance is

```math
R_\phi = \rho_{20}\left[1+\alpha_{20}(T-20^\circ\mathrm C)\right]
\frac{\ell_\phi}{A_\mathrm{Cu}} + R_\mathrm{joint}.
```

The default temperature coefficient, 0.00393/K, is the International Annealed
Copper Standard value near 20 C. Insulated wire diameter controls packing;
bare AWG area controls resistance. Leads, joints, and parallel paths remain
separate inputs.

For `inductance_model=:coreless_filament`, every nested round-wire turn in all
series phase coils is represented by straight line elements. The phase partial
inductance follows the Neumann integral

```math
L_\phi = \frac{\mu_0}{4\pi}\sum_{i,j}
\oint_{C_i}\oint_{C_j}\frac{d\boldsymbol\ell_i\cdot d\boldsymbol\ell_j}
{|\mathbf r_i-\mathbf r_j|}.
```

The round-wire self term uses the geometric-mean-distance radius
`r*exp(-1/4)`. This is a free-space conductor calculation, so magnetic backing
and independently justified lead inductance must be supplied separately. The
path subdivision should be converged for the geometry; it is not a calibration
parameter. Slotless-winding inductance is particularly sensitive to coil
geometry and end/turn leakage, as discussed by Jumayev et al.,
<https://doi.org/10.1108/COMPEL-08-2014-0207>.

For `halbach_field_model=:segmented_cuboid`, each finite, uniformly magnetized
block is evaluated from the closed-form field of its rectangular bound
surface charges. The model sums all blocks, numerically integrates normal
field over each packed turn aperture as the rotor is advanced, and
Fourier-projects the resulting phase linkage. This avoids a
single mean-radius `B*A*N` approximation and represents actual segmentation,
block dimensions, remanence, magnetization order, coil shape, axial layers,
and magnet-to-wire distance. The surface-charge treatment is consistent with
polyhedral permanent-magnet field formulations such as O'Connell,
<https://doi.org/10.1016/j.jmmm.2020.166894>, while spatial coil integration is
the essential step in axial-flux EMF prediction (for example,
<https://www.mdpi.com/2079-9292/14/14/2901>).

The current cuboid implementation is deliberately a single-rotor free-space
model: it assumes unity recoil permeability and omits magnetic backing,
demagnetization, adhesive gaps, placement/orientation tolerances, and eddy
currents. Use a validated subdomain solution or 3-D FEA when those effects are
material. A geometry-only result should be validated with open-circuit EMF,
four-wire resistance, and an independent phase-inductance measurement before
being treated as predictive design evidence.
