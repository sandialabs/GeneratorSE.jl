# Examples

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

Run the axial Halbach validation example:

```bash
julia --project examples/pmsg_axial_halbach_validation.jl
```

This saves:

- `example_pmsg_axial_halbach_flux_validation.pdf`
- `example_pmsg_axial_halbach_generator_validation.pdf`
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
