# Examples

The axial and arms examples use GeneratorSE's corrected electrical conventions:
`p` is pole pairs, so `omega_e = p * omega_m`, and the current estimate divides
the reactive term by `omega_e * L_s`. This differs from the original translated
code and is documented in the package README and theory docs.

Run the axial operating-surface example from the repository root:

```bash
julia --project examples/pmsg_axial_operating_surfaces.jl
```

Generated figures are saved to `figs/` as transparent PDFs:

- `example_pmsg_axial_efficiency_surface.pdf`
- `example_pmsg_axial_losses_surface.pdf`
- `example_pmsg_axial_current_surface.pdf`
- `example_pmsg_axial_voltage_surface.pdf`

Optional controls:

```bash
GENERATORSE_EXAMPLE_NRPM=25 GENERATORSE_EXAMPLE_NTORQUE=25 julia --project examples/pmsg_axial_operating_surfaces.jl
```

Run the axial Halbach regression/sanity example:

```bash
julia --project examples/pmsg_axial_halbach_validation.jl
```

This is not an independent validation set: the ideal Halbach target mirrors the
implementation equation, and the other comparisons are normalized sanity checks.
This saves:

- `example_pmsg_axial_halbach_flux_regression.pdf`
- `example_pmsg_axial_halbach_generator_sanity.pdf`
- `example_pmsg_axial_halbach_metrics.csv`

For hand-wound or coreless prototypes, the Halbach model accepts physical
winding inputs without requiring measured electrical constants. Resistance uses
the following precedence:

1. `phase_resistance`;
2. any available physical inputs: `conductor_area` (or `wire_gauge_awg`) and
   `mean_turn_length`, or the complete trapezoid support geometry
   `coil_inner_radius`, `coil_outer_radius`, and `coil_span_angle`; and
3. the legacy slot-fill estimate when none of those inputs is supplied.

`phase_lead_length`, `parallel_paths`, and the existing `resist_Cu` input make
the remaining resistance assumptions explicit. `mean_turn_length` is the full
length of one turn, not one coil side.
Wire area and turn geometry may be supplied independently; an omitted quantity
uses the model's internal winding-geometry or slot-fill assumption.

Inductance similarly uses `phase_inductance` first. Otherwise,
`inductance_model=:concentrated_coreless` requires `turns_per_coil`,
`coils_in_series_per_phase`, and `inductance_coil_area`; optional mutual and
leakage assumptions are supplied with `coil_mutual_coupling` and
`phase_leakage_inductance`. Omitting these inputs preserves the legacy
distributed-winding estimate.

Flux linkage can use an explicit `phase_flux_linkage`, or an
`airgap_flux_density` with `effective_flux_area`, `winding_factor`, and
`flux_linkage_factor`. Here `effective_flux_area` means the field-weighted area
linked by one coil/pole, not necessarily the gross coil envelope. If no explicit
flux input is given, the Halbach field and annular-area defaults are retained.

Run the ORNL Prius PMSM performance-map validation example:

```bash
julia --project examples/pmsm_prius_performance_map_validation.jl
```

This uses embedded, documented points digitized from the published ORNL 2010
Prius MG2 PMSM speed-torque efficiency maps and torque-current data.

This saves:

- `example_pmsm_prius_performance_map_efficiency.pdf`
- `example_pmsm_prius_performance_map_electrical.pdf`
- `example_pmsm_prius_performance_map_line_comparison.pdf`
- `example_pmsm_prius_performance_map_points.csv`
- `example_pmsm_prius_performance_map_metrics.csv`

Run a realistic dynamic-response example:

```bash
julia --project examples/pmsg_axial_dynamic_response.jl
```

This saves:

- `example_pmsg_axial_dynamic_response.pdf`
