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

Run the axial Halbach final-quantity validation example:

```bash
julia --project examples/pmsg_axial_halbach_final_quantity_validation.jl
```

This calls `PMSG_axial_Halbach` on measured nonzero-speed axial-flux PM motor
operating points and compares torque/speed against current and voltage.

This saves:

- `example_pmsg_axial_halbach_final_quantity_validation.pdf`
- `example_pmsg_axial_halbach_final_quantity_points.csv`
- `example_pmsg_axial_halbach_final_quantity_metrics.csv`

Run the Sandia WEC Spar V0 prototype Halbach model:

```bash
julia --project examples/pmsg_axial_halbach_sandia_wec_prototype.jl
```

This uses the local milestone DOCX and `Context.md` assumptions to model the
small coreless prototype with explicit turns, coil area, and measured phase R/L.

This saves:

- `sandia_wec_prototype_model.csv`

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
