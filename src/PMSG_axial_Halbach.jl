# Axial-flux variant with a screening-level Halbach rotor model.
# The output tuple mirrors PMSG_axial so downstream sizing code can swap models.

function _sinc_unity(x)
    if abs(x) < 1.0e-6
        x2 = x * x
        return one(x) - x2 / 6 + x2 * x2 / 120
    end
    return sin(x) / x
end

function _awg_conductor_area(awg)
    diameter_m = 0.005 * 0.0254 * 92^((36 - awg) / 39)
    return pi * diameter_m^2 / 4
end

"""
    copper_resistivity_at_temperature(rho_20, temperature_c; alpha_20=0.00393)

Return copper resistivity at `temperature_c` from a 20 °C reference.  The
default temperature coefficient is the International Annealed Copper Standard
value for annealed winding wire near 20 °C.
"""
function copper_resistivity_at_temperature(rho_20, temperature_c; alpha_20 = 0.00393)
    return rho_20 * (1 + alpha_20 * (temperature_c - 20))
end

function _gauss_legendre(order)
    if order == 2
        a = inv(sqrt(3.0))
        return (-a, a), (1.0, 1.0)
    elseif order == 3
        a = sqrt(3 / 5)
        return (-a, 0.0, a), (5 / 9, 8 / 9, 5 / 9)
    elseif order == 4
        return (
            -0.8611363115940526,
            -0.3399810435848563,
            0.3399810435848563,
            0.8611363115940526,
        ), (
            0.3478548451374538,
            0.6521451548625461,
            0.6521451548625461,
            0.3478548451374538,
        )
    elseif order == 6
        return (
            -0.9324695142031521,
            -0.6612093864662645,
            -0.2386191860831969,
            0.2386191860831969,
            0.6612093864662645,
            0.9324695142031521,
        ), (
            0.1713244923791704,
            0.3607615730481386,
            0.4679139345726910,
            0.4679139345726910,
            0.3607615730481386,
            0.1713244923791704,
        )
    end
    throw(ArgumentError("quadrature order must be 2, 3, 4, or 6"))
end

_vadd(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
_vscale(a, s) = (a[1] * s, a[2] * s, a[3] * s)
_vdot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

_vcross(a, b) = (
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1],
)

function _append_magnetic_charge_face!(sources, center, axis_u, axis_v, length_u, length_v, charge_fraction, B_r)
    abs(charge_fraction) <= 1.0e-14 && return sources
    push!(sources, (
        center = center,
        axis_u = axis_u,
        axis_v = axis_v,
        half_u = length_u / 2,
        half_v = length_v / 2,
        coefficient = B_r * charge_fraction / (4 * pi),
    ))
    return sources
end

function _rectangular_charge_face_Bz(point, face)
    # Closed-form electric/gravitational field integral for a uniformly
    # charged rectangle, here multiplied by Br/(4*pi) to obtain B.  The local
    # coordinates are right-handed (u, v, n=u×v).  Using asinh for the
    # in-plane components avoids cancellation in log(v + R) near an edge.
    axis_n = _vcross(face.axis_u, face.axis_v)
    offset = (
        point[1] - face.center[1],
        point[2] - face.center[2],
        point[3] - face.center[3],
    )
    x = _vdot(offset, face.axis_u)
    y = _vdot(offset, face.axis_v)
    z = _vdot(offset, axis_n)
    us = (x + face.half_u, x - face.half_u)
    vs = (y + face.half_v, y - face.half_v)
    signs = (one(x), -one(x))
    integral_u = zero(x)
    integral_v = zero(x)
    integral_n = zero(x)
    for i in 1:2, j in 1:2
        u = us[i]
        v = vs[j]
        sign = signs[i] * signs[j]
        radius = sqrt(u * u + v * v + z * z)
        integral_u -= sign * asinh(v / sqrt(u * u + z * z))
        integral_v -= sign * asinh(u / sqrt(v * v + z * z))
        if abs(axis_n[3]) > 1.0e-14
            integral_n += sign * atan(u * v / (z * radius))
        end
    end
    return face.coefficient * (
        integral_u * face.axis_u[3] +
        integral_v * face.axis_v[3] +
        integral_n * axis_n[3]
    )
end

"""
    segmented_halbach_winding_properties(; kwargs...)

Calculate the free-space field and phase fundamental flux linkage of a
single-rotor segmented axial Halbach array.  Each uniformly magnetized cuboid
is represented by the closed-form field of its rectangular bound magnetic
surface charges.
The field is integrated over the actual annular-sector coil apertures at a set
of rotor positions; the phase fundamental is then obtained by Fourier
projection.  This replaces the lumped `B*A*N*k_w` approximation when selected
from `PMSG_axial_Halbach` with `halbach_field_model=:segmented_cuboid`.
When `wire_outer_diameter` and `turns_per_layer` are supplied, every nested
turn aperture and axial layer is integrated separately; otherwise the supplied
coil aperture is multiplied by `turns_per_coil`.
For `winding_geometry_reference=:inner_support`, conductor paths use exact
parallel offsets while the field integral uses an equal-area annular-sector
proxy for each rounded offset aperture.

The standard sequence is `+z, -tangential, -z, +tangential`, matching a
four-step Halbach array whose strong side is above the rotor.  Back iron and
non-unity recoil permeability are not included; use 3-D FEA or an appropriate
subdomain model when those effects are material.
"""
function segmented_halbach_winding_properties(;
    pole_pairs,
    magnet_inner_radius,
    magnet_outer_radius,
    magnet_thickness,
    magnet_tangential_width,
    magnets_total,
    B_r,
    winding_plane_offset,
    coil_inner_radius,
    coil_outer_radius,
    coil_span_angle,
    turns_per_coil,
    coils_in_series_per_phase,
    phase_coil_offset = 0.0,
    wire_outer_diameter = nothing,
    turns_per_layer = nothing,
    turn_pitch = nothing,
    layer_pitch = nothing,
    winding_geometry_reference = :mean_turn,
    winding_support_clearance = 0.0,
    phase_count = 3,
    magnetization_rotation = -1.0,
    coil_quadrature_order = 6,
    rotor_samples = 48,
)
    p = pole_pairs
    n_magnets = Int(round(magnets_total))
    n_coils = Int(round(coils_in_series_per_phase))
    n_samples = Int(round(rotor_samples))
    if n_magnets < 4 || n_coils < 1 || n_samples < 8
        throw(ArgumentError("segmented Halbach model requires at least 4 magnets, 1 phase coil, and 8 rotor samples"))
    end
    if magnet_outer_radius <= magnet_inner_radius || magnet_thickness <= 0 || magnet_tangential_width <= 0
        throw(ArgumentError("magnet dimensions must be positive and ordered"))
    end
    if coil_outer_radius <= coil_inner_radius || coil_inner_radius < 0 || coil_span_angle <= 0
        throw(ArgumentError("coil annular-sector geometry must be positive and ordered"))
    end
    if winding_plane_offset <= 0 || turns_per_coil <= 0
        throw(ArgumentError("winding-plane offset and turns per coil must be positive"))
    end

    T = promote_type(
        typeof(p), typeof(magnet_inner_radius), typeof(magnet_outer_radius),
        typeof(magnet_thickness), typeof(magnet_tangential_width), typeof(B_r),
        typeof(winding_plane_offset), typeof(coil_inner_radius),
        typeof(coil_outer_radius), typeof(coil_span_angle), typeof(turns_per_coil),
        typeof(phase_coil_offset), typeof(magnetization_rotation), Float64,
    )
    Face = NamedTuple{
        (:center, :axis_u, :axis_v, :half_u, :half_v, :coefficient),
        Tuple{NTuple{3,T},NTuple{3,T},NTuple{3,T},T,T,T},
    }
    sources = Face[]
    radius = convert(T, 0.5) * (magnet_inner_radius + magnet_outer_radius)
    radial_length = magnet_outer_radius - magnet_inner_radius
    z_center = -convert(T, 0.5) * magnet_thickness
    two_pi = convert(T, 2 * pi)

    for j = 0:(n_magnets - 1)
        theta = two_pi * convert(T, j) / convert(T, n_magnets)
        er = (cos(theta), sin(theta), zero(T))
        et = (-sin(theta), cos(theta), zero(T))
        ez = (zero(T), zero(T), one(T))
        center = (radius * er[1], radius * er[2], z_center)

        magnetization_phase = two_pi * p * convert(T, j) / convert(T, n_magnets)
        m_t = magnetization_rotation * sin(magnetization_phase)
        m_z = cos(magnetization_phase)

        for face_sign in (-1, 1)
            s = convert(T, face_sign)
            tangential_center = _vadd(center, _vscale(et, s * magnet_tangential_width / 2))
            _append_magnetic_charge_face!(
                sources, tangential_center, er, ez, radial_length, magnet_thickness,
                s * m_t, B_r,
            )
            axial_center = _vadd(center, _vscale(ez, s * magnet_thickness / 2))
            _append_magnetic_charge_face!(
                sources, axial_center, er, et, radial_length, magnet_tangential_width,
                s * m_z, B_r,
            )
        end
    end

    function Bz_at(r, theta)
        x = r * cos(theta)
        y = r * sin(theta)
        z = winding_plane_offset
        Bz = zero(T)
        for source in sources
            Bz += _rectangular_charge_face_Bz((x, y, z), source)
        end
        return Bz
    end

    packed_turns = if wire_outer_diameter !== nothing || turns_per_layer !== nothing
        if wire_outer_diameter === nothing || turns_per_layer === nothing
            throw(ArgumentError("wire_outer_diameter and turns_per_layer must be provided together"))
        end
        _wedge_turn_geometries(
            coil_inner_radius,
            coil_outer_radius,
            coil_span_angle,
            turns_per_coil,
            turns_per_layer,
            turn_pitch === nothing ? wire_outer_diameter : turn_pitch,
            layer_pitch === nothing ? wire_outer_diameter : layer_pitch,
            winding_geometry_reference,
            winding_support_clearance,
        )
    else
        nothing
    end
    winding_clearance = packed_turns === nothing ? nothing :
        _validate_interleaved_winding_clearance(
            packed_turns, phase_count, coils_in_series_per_phase,
            wire_outer_diameter,
        )

    coil_nodes, coil_weights = _gauss_legendre(coil_quadrature_order)
    function coil_flux(center_angle, r_inner, r_outer, span, z_offset)
        flux = zero(T)
        radial_jacobian = (r_outer - r_inner) / 2
        angular_jacobian = span / 2
        radial_center = (r_outer + r_inner) / 2
        for (xr, wr) in zip(coil_nodes, coil_weights), (xt, wt) in zip(coil_nodes, coil_weights)
            r = radial_center + radial_jacobian * xr
            theta = center_angle + angular_jacobian * xt
            # Axial layers see slightly different field.  Re-evaluate the
            # magnetic surface-charge sum at the center plane of each turn.
            x = r * cos(theta)
            y = r * sin(theta)
            z = winding_plane_offset + z_offset
            Bz = zero(T)
            for source in sources
                Bz += _rectangular_charge_face_Bz((x, y, z), source)
            end
            flux += Bz * r * radial_jacobian * angular_jacobian * wr * wt
        end
        return flux
    end

    rotor_angles = Vector{T}(undef, n_samples)
    linkage = Vector{T}(undef, n_samples)
    field = Vector{T}(undef, n_samples)
    coil_pitch = two_pi / convert(T, n_coils)
    for sample = 1:n_samples
        rotor_angle = two_pi * convert(T, sample - 1) / (p * convert(T, n_samples))
        rotor_angles[sample] = rotor_angle
        phase_flux = zero(T)
        for coil = 0:(n_coils - 1)
            coil_angle = phase_coil_offset + convert(T, coil) * coil_pitch - rotor_angle
            if packed_turns === nothing
                phase_flux += turns_per_coil * coil_flux(
                    coil_angle,
                    coil_inner_radius,
                    coil_outer_radius,
                    coil_span_angle,
                    zero(T),
                )
            else
                for turn in packed_turns
                    phase_flux += coil_flux(coil_angle, turn.r_in, turn.r_out, turn.span, turn.z)
                end
            end
        end
        linkage[sample] = phase_flux
        field[sample] = Bz_at(radius, phase_coil_offset - rotor_angle)
    end

    lambda_cos = zero(T)
    lambda_sin = zero(T)
    field_cos = zero(T)
    field_sin = zero(T)
    for sample = 1:n_samples
        electrical_angle = p * rotor_angles[sample]
        c = cos(electrical_angle)
        s = sin(electrical_angle)
        lambda_cos += linkage[sample] * c
        lambda_sin += linkage[sample] * s
        field_cos += field[sample] * c
        field_sin += field[sample] * s
    end
    scale = convert(T, 2) / convert(T, n_samples)
    lambda_fundamental = hypot(scale * lambda_cos, scale * lambda_sin)
    field_fundamental = hypot(scale * field_cos, scale * field_sin)

    return (
        B_fundamental = field_fundamental,
        B_peak_at_mean_radius = maximum(abs, field),
        phase_flux_linkage_fundamental = lambda_fundamental,
        rotor_angles = rotor_angles,
        phase_flux_linkage = linkage,
        mean_radius_field = field,
        magnetic_charge_sources = length(sources),
        packed_turn_geometry = packed_turns !== nothing,
        winding_clearance = winding_clearance,
    )
end

function _wedge_turn_geometries(
    coil_inner_radius,
    coil_outer_radius,
    coil_span_angle,
    turns_per_coil,
    turns_per_layer,
    turn_pitch,
    layer_pitch,
    geometry_reference,
    support_clearance = 0.0,
)
    n_turns = Int(round(turns_per_coil))
    n_per_layer = Int(round(turns_per_layer))
    if n_turns < 1 || n_per_layer < 1 || turn_pitch <= 0 || layer_pitch <= 0 || support_clearance < 0
        throw(ArgumentError("turn count, turns per layer, and winding pitches must be positive"))
    end
    n_layers = cld(n_turns, n_per_layer)
    radius_mean = (coil_inner_radius + coil_outer_radius) / 2
    T = promote_type(
        typeof(coil_inner_radius), typeof(coil_outer_radius), typeof(coil_span_angle),
        typeof(turn_pitch), typeof(layer_pitch), Float64,
    )
    geometries = NamedTuple[]
    for turn = 0:(n_turns - 1)
        layer = fld(turn, n_per_layer)
        slot = mod(turn, n_per_layer)
        turns_this_layer = min(n_per_layer, n_turns - layer * n_per_layer)
        inplane_offset = if geometry_reference == :mean_turn
            (convert(T, slot) - convert(T, turns_this_layer - 1) / 2) * turn_pitch
        elseif geometry_reference == :inner_support
            support_clearance + (convert(T, slot) + convert(T, 0.5)) * turn_pitch
        else
            throw(ArgumentError("winding_geometry_reference must be :mean_turn or :inner_support"))
        end
        axial_offset = (convert(T, layer) - convert(T, n_layers - 1) / 2) * layer_pitch
        r_in = coil_inner_radius - inplane_offset
        r_out = coil_outer_radius + inplane_offset
        # For an inner-support reference, preserve the exact area of a
        # parallel offset of the bobbin aperture (Steiner formula) in the
        # annular-sector proxy used by the magnetic-field quadrature. The
        # inductance path itself is generated as the exact offset curve below.
        span = if geometry_reference == :inner_support
            support_area = 0.5 * (coil_outer_radius^2 - coil_inner_radius^2) * coil_span_angle
            support_perimeter = 2 * (coil_outer_radius - coil_inner_radius) +
                                coil_span_angle * (coil_inner_radius + coil_outer_radius)
            offset_area = support_area + support_perimeter * inplane_offset + pi * inplane_offset^2
            2 * offset_area / (r_out^2 - r_in^2)
        else
            coil_span_angle + 2 * inplane_offset / radius_mean
        end
        if r_in <= 0 || r_out <= r_in || span <= 0
            throw(ArgumentError("wire packing creates an invalid turn geometry"))
        end
        push!(geometries, (
            r_in = r_in,
            r_out = r_out,
            span = span,
            z = axial_offset,
            support_r_in = coil_inner_radius,
            support_r_out = coil_outer_radius,
            support_span = coil_span_angle,
            inplane_offset = inplane_offset,
            exact_parallel_offset = geometry_reference == :inner_support,
        ))
    end
    return geometries
end

function _wedge_turn_length(turn)
    if turn.exact_parallel_offset
        support_perimeter = 2 * (turn.support_r_out - turn.support_r_in) +
                            turn.support_span * (turn.support_r_in + turn.support_r_out)
        return support_perimeter + 2 * pi * turn.inplane_offset
    end
    return 2 * (turn.r_out - turn.r_in) + turn.span * (turn.r_in + turn.r_out)
end

function _validate_interleaved_winding_clearance(
    turns, phase_count, coils_per_phase, wire_outer_diameter,
)
    isempty(turns) && return nothing
    if !turns[1].exact_parallel_offset
        return nothing
    end
    n_phases = Int(round(phase_count))
    n_coils = Int(round(coils_per_phase))
    if n_phases < 2 || abs(phase_count - n_phases) > 1.0e-8
        throw(ArgumentError("phase_count must be an integer of at least two"))
    end
    if n_coils < 1 || abs(coils_per_phase - n_coils) > 1.0e-8
        throw(ArgumentError("coils_per_phase must be a positive integer"))
    end
    phase_pitch = 2 * pi / (n_phases * n_coils)
    maximum_center_offset = maximum(turn.inplane_offset for turn in turns)
    outer_winding_offset = maximum_center_offset + wire_outer_diameter / 2
    occupied_span = turns[1].support_span +
                    2 * atan(outer_winding_offset / turns[1].support_r_in)
    margin = phase_pitch - occupied_span
    if margin <= 0
        throw(ArgumentError(
            "inner-support turn packing overlaps the adjacent phase coil; " *
            "reduce turns_per_layer, wire diameter/clearance, or coil support span",
        ))
    end
    return (
        phase_coil_pitch_angle = phase_pitch,
        occupied_coil_span_angle = occupied_span,
        coil_packing_margin_angle = margin,
    )
end

function _append_line_elements!(elements, p0, p1, conductor, subdivisions)
    n = Int(round(subdivisions))
    for k = 0:(n - 1)
        a = k / n
        b = (k + 1) / n
        start = (
            p0[1] + a * (p1[1] - p0[1]),
            p0[2] + a * (p1[2] - p0[2]),
            p0[3] + a * (p1[3] - p0[3]),
        )
        stop = (
            p0[1] + b * (p1[1] - p0[1]),
            p0[2] + b * (p1[2] - p0[2]),
            p0[3] + b * (p1[3] - p0[3]),
        )
        dl = (stop[1] - start[1], stop[2] - start[2], stop[3] - start[3])
        midpoint = ((start[1] + stop[1]) / 2, (start[2] + stop[2]) / 2, (start[3] + stop[3]) / 2)
        push!(elements, (r = midpoint, dl = dl, conductor = conductor))
    end
    return elements
end

function _append_arc_elements!(elements, radius, theta0, theta1, z, conductor, subdivisions)
    n = Int(round(subdivisions))
    for k = 0:(n - 1)
        a0 = theta0 + (theta1 - theta0) * k / n
        a1 = theta0 + (theta1 - theta0) * (k + 1) / n
        p0 = (radius * cos(a0), radius * sin(a0), z)
        p1 = (radius * cos(a1), radius * sin(a1), z)
        _append_line_elements!(elements, p0, p1, conductor, 1)
    end
    return elements
end

function _append_centered_arc_elements!(
    elements, center, radius, theta0, theta1, z, conductor, subdivisions,
)
    radius <= 0 && return elements
    n = Int(round(subdivisions))
    for k = 0:(n - 1)
        a0 = theta0 + (theta1 - theta0) * k / n
        a1 = theta0 + (theta1 - theta0) * (k + 1) / n
        p0 = (center[1] + radius * cos(a0), center[2] + radius * sin(a0), z)
        p1 = (center[1] + radius * cos(a1), center[2] + radius * sin(a1), z)
        _append_line_elements!(elements, p0, p1, conductor, 1)
    end
    return elements
end

function _append_parallel_offset_wedge_loop_elements!(
    elements, turn, center_angle, conductor, subdivisions,
)
    half_span = turn.support_span / 2
    theta_low = center_angle - half_span
    theta_high = center_angle + half_span
    r_in = turn.support_r_in
    r_out = turn.support_r_out
    offset = turn.inplane_offset
    z = turn.z
    er_low = (cos(theta_low), sin(theta_low), zero(z))
    et_low = (-sin(theta_low), cos(theta_low), zero(z))
    er_high = (cos(theta_high), sin(theta_high), zero(z))
    et_high = (-sin(theta_high), cos(theta_high), zero(z))
    inner_low = (r_in * er_low[1], r_in * er_low[2])
    outer_low = (r_out * er_low[1], r_out * er_low[2])
    inner_high = (r_in * er_high[1], r_in * er_high[2])
    outer_high = (r_out * er_high[1], r_out * er_high[2])

    p_inner_low = (
        inner_low[1] - offset * et_low[1],
        inner_low[2] - offset * et_low[2],
        z,
    )
    p_outer_low = (
        outer_low[1] - offset * et_low[1],
        outer_low[2] - offset * et_low[2],
        z,
    )
    p_outer_high = (
        outer_high[1] + offset * et_high[1],
        outer_high[2] + offset * et_high[2],
        z,
    )
    p_inner_high = (
        inner_high[1] + offset * et_high[1],
        inner_high[2] + offset * et_high[2],
        z,
    )

    _append_line_elements!(elements, p_inner_low, p_outer_low, conductor, subdivisions)
    _append_centered_arc_elements!(
        elements, outer_low, offset, theta_low - pi / 2, theta_low,
        z, conductor, subdivisions,
    )
    _append_arc_elements!(elements, r_out + offset, theta_low, theta_high, z, conductor, subdivisions)
    _append_centered_arc_elements!(
        elements, outer_high, offset, theta_high, theta_high + pi / 2,
        z, conductor, subdivisions,
    )
    _append_line_elements!(elements, p_outer_high, p_inner_high, conductor, subdivisions)
    _append_centered_arc_elements!(
        elements, inner_high, offset, theta_high + pi / 2, theta_high + pi,
        z, conductor, subdivisions,
    )
    _append_arc_elements!(elements, r_in - offset, theta_high, theta_low, z, conductor, subdivisions)
    _append_centered_arc_elements!(
        elements, inner_low, offset, theta_low + pi, theta_low + 3 * pi / 2,
        z, conductor, subdivisions,
    )
    return elements
end

function _append_wedge_loop_elements!(elements, turn, center_angle, conductor, subdivisions)
    if turn.exact_parallel_offset
        return _append_parallel_offset_wedge_loop_elements!(
            elements, turn, center_angle, conductor, subdivisions,
        )
    end
    half_span = turn.span / 2
    theta_low = center_angle - half_span
    theta_high = center_angle + half_span
    p_inner_low = (turn.r_in * cos(theta_low), turn.r_in * sin(theta_low), turn.z)
    p_outer_low = (turn.r_out * cos(theta_low), turn.r_out * sin(theta_low), turn.z)
    p_outer_high = (turn.r_out * cos(theta_high), turn.r_out * sin(theta_high), turn.z)
    p_inner_high = (turn.r_in * cos(theta_high), turn.r_in * sin(theta_high), turn.z)
    _append_line_elements!(elements, p_inner_low, p_outer_low, conductor, subdivisions)
    _append_arc_elements!(elements, turn.r_out, theta_low, theta_high, turn.z, conductor, subdivisions)
    _append_line_elements!(elements, p_outer_high, p_inner_high, conductor, subdivisions)
    _append_arc_elements!(elements, turn.r_in, theta_high, theta_low, turn.z, conductor, subdivisions)
    return elements
end

function _coil_winding_elements(turns, center_angle, path_subdivisions, ::Type{T}) where {T}
    Element = NamedTuple{(:r, :dl, :conductor),Tuple{NTuple{3,T},NTuple{3,T},Int}}
    elements = Element[]
    for (conductor, turn) in enumerate(turns)
        _append_wedge_loop_elements!(
            elements, turn, center_angle, conductor, path_subdivisions,
        )
    end
    return elements
end

function _coil_self_inductance(elements, gmd_radius, coefficient)
    L = zero(coefficient)
    @inbounds for i in eachindex(elements)
        ei = elements[i]
        length_dl = sqrt(_vdot(ei.dl, ei.dl))
        self_double_integral = 2 * (
            length_dl * asinh(length_dl / gmd_radius) -
            sqrt(length_dl^2 + gmd_radius^2) + gmd_radius
        )
        L += coefficient * self_double_integral
        for j = (i + 1):length(elements)
            ej = elements[j]
            dx = ei.r[1] - ej.r[1]
            dy = ei.r[2] - ej.r[2]
            dz = ei.r[3] - ej.r[3]
            distance2 = dx * dx + dy * dy + dz * dz
            # Adjacent pieces of one physical conductor share the same
            # round-wire GMD regularization as the analytic self term.
            distance = ei.conductor == ej.conductor ?
                sqrt(distance2 + gmd_radius^2) : sqrt(distance2)
            L += 2 * coefficient * _vdot(ei.dl, ej.dl) / distance
        end
    end
    return L
end

function _coil_mutual_inductance(elements_a, elements_b, coefficient)
    M = zero(coefficient)
    @inbounds for ea in elements_a, eb in elements_b
        dx = ea.r[1] - eb.r[1]
        dy = ea.r[2] - eb.r[2]
        dz = ea.r[3] - eb.r[3]
        M += coefficient * _vdot(ea.dl, eb.dl) / sqrt(dx * dx + dy * dy + dz * dz)
    end
    return M
end

"""
    coreless_winding_inductance(; kwargs...)

Calculate the positive-sequence (controller `dq`) inductance from the Neumann
partial-inductance integral over the actual series-connected air-core turn
paths. Individual nested turns, axial layers, finite round-wire self terms,
same-phase coupling, and cross-phase mutual coupling are included. Phases are
assumed to be uniformly interleaved around the stator. The returned
`phase_inductance_matrix` is the full phase-variable matrix; for a balanced
three-phase winding, `phase_inductance == dq_inductance == L_aa - M_ab` and a
two-terminal line-to-line test gives `line_line_inductance == 2*dq_inductance`.

The result excludes the phase jumpers and magnetic backing. An independently
calculated, uncoupled per-phase lead contribution can be added with
`phase_lead_inductance`. Permanent-magnet remanence is deliberately absent:
it sets source flux linkage, not the incremental inductance of an ironless
winding.
"""
function coreless_winding_inductance(;
    coil_inner_radius,
    coil_outer_radius,
    coil_span_angle,
    turns_per_coil,
    coils_in_series_per_phase,
    wire_diameter,
    conductor_diameter = wire_diameter,
    turns_per_layer,
    turn_pitch = wire_diameter,
    layer_pitch = wire_diameter,
    winding_geometry_reference = :mean_turn,
    winding_support_clearance = 0.0,
    path_subdivisions = 12,
    phase_count = 3,
    phase_lead_inductance = 0.0,
    mu_0 = 4 * pi * 1.0e-7,
)
    if wire_diameter <= 0 || conductor_diameter <= 0 || conductor_diameter > wire_diameter || phase_lead_inductance < 0
        throw(ArgumentError("wire diameters must be positive, bare conductor no larger than the insulated wire, and lead inductance nonnegative"))
    end
    turns = _wedge_turn_geometries(
        coil_inner_radius, coil_outer_radius, coil_span_angle, turns_per_coil,
        turns_per_layer, turn_pitch, layer_pitch, winding_geometry_reference,
        winding_support_clearance,
    )
    n_coils = Int(round(coils_in_series_per_phase))
    n_coils < 1 && throw(ArgumentError("coils in series per phase must be positive"))
    n_phases = Int(round(phase_count))
    if n_phases < 2 || abs(phase_count - n_phases) > 1.0e-8
        throw(ArgumentError("phase_count must be an integer of at least two"))
    end
    winding_clearance = _validate_interleaved_winding_clearance(
        turns, n_phases, n_coils, wire_diameter,
    )
    T = promote_type(
        typeof(coil_inner_radius), typeof(coil_outer_radius), typeof(coil_span_angle),
        typeof(wire_diameter), typeof(conductor_diameter), typeof(mu_0), Float64,
    )
    gmd_radius = 0.5 * conductor_diameter * exp(-0.25)
    coefficient = mu_0 / (4 * pi)
    reference_elements = _coil_winding_elements(turns, zero(T), path_subdivisions, T)
    coil_self = _coil_self_inductance(reference_elements, gmd_radius, coefficient)

    # Rotational symmetry avoids constructing every turn in every phase while
    # retaining the complete phase matrix. For one phase, n_coils times the
    # coupling from a reference coil to all other same-phase coils is exactly
    # the ordered-pair sum that appears in magnetic co-energy.
    phase_self_winding = n_coils * coil_self
    for coil = 1:(n_coils - 1)
        elements = _coil_winding_elements(
            turns, 2 * pi * coil / n_coils, path_subdivisions, T,
        )
        phase_self_winding += n_coils * _coil_mutual_inductance(
            reference_elements, elements, coefficient,
        )
    end

    # A uniformly interleaved m-phase winding has phase-k coils shifted by
    # k/(m*n_coils) mechanical revolutions from phase 1. Each row of L_phase
    # is circulant, so only the reference-phase mutuals are required.
    phase_mutual = Vector{T}(undef, n_phases - 1)
    for phase = 1:(n_phases - 1)
        M = zero(T)
        phase_offset = 2 * pi * phase / (n_phases * n_coils)
        for coil = 0:(n_coils - 1)
            elements = _coil_winding_elements(
                turns, phase_offset + 2 * pi * coil / n_coils,
                path_subdivisions, T,
            )
            M += n_coils * _coil_mutual_inductance(
                reference_elements, elements, coefficient,
            )
        end
        phase_mutual[phase] = M
    end
    for phase = 1:(n_phases - 1)
        conjugate_phase = n_phases - phase
        mutual_average = (phase_mutual[phase] + phase_mutual[conjugate_phase]) / 2
        phase_mutual[phase] = mutual_average
        phase_mutual[conjugate_phase] = mutual_average
    end

    phase_self = phase_self_winding + phase_lead_inductance
    phase_inductance_matrix = Matrix{T}(undef, n_phases, n_phases)
    for row = 1:n_phases, column = 1:n_phases
        separation = mod(column - row, n_phases)
        phase_inductance_matrix[row, column] = separation == 0 ?
            phase_self : phase_mutual[separation]
    end

    # The positive-sequence eigenvalue is the stationary alpha-beta / rotating
    # dq inductance. Reciprocity makes the imaginary Fourier component vanish;
    # averaging conjugate mutual terms also suppresses quadrature roundoff.
    dq_winding_inductance = phase_self_winding
    for phase = 1:(n_phases - 1)
        dq_winding_inductance += phase_mutual[phase] * cos(2 * pi * phase / n_phases)
    end
    dq_inductance = dq_winding_inductance + phase_lead_inductance
    line_line_inductance = n_phases == 3 ? 2 * dq_inductance : nothing
    phase_turn_length = n_coils * sum(_wedge_turn_length(turn) for turn in turns)
    return (
        phase_inductance = dq_inductance,
        winding_inductance = dq_winding_inductance,
        dq_inductance = dq_inductance,
        phase_self_inductance = phase_self,
        phase_mutual_inductances = phase_mutual,
        phase_inductance_matrix = phase_inductance_matrix,
        line_line_inductance = line_line_inductance,
        phase_lead_inductance = phase_lead_inductance,
        conductor_diameter = conductor_diameter,
        wire_outer_diameter = wire_diameter,
        phase_turn_length = phase_turn_length,
        phase_count = n_phases,
        winding_clearance = winding_clearance,
        turns_per_layer = turns_per_layer,
        axial_layers = cld(Int(round(turns_per_coil)), Int(round(turns_per_layer))),
        conductor_loops = n_coils * length(turns),
        line_elements = n_coils * length(reference_elements),
    )
end

function _halbach_segmentation_factor(segment_count, magnet_width_ratio, field_model)
    x_segment = pi / (2 * segment_count)

    if field_model == :ideal_sheet
        return _sinc_unity(x_segment)
    elseif field_model == :finite_width_harmonic
        segment_fill = _smooth_min(_smooth_max(magnet_width_ratio, 0.0), 1.0)
        return segment_fill * _sinc_unity(segment_fill * x_segment)
    end

    throw(ArgumentError("halbach_field_model must be :ideal_sheet or :finite_width_harmonic"))
end

function halbach_fundamental_flux_density(
    B_r,
    h_m,
    len_ag,
    tau_p;
    mu_r = 1.06,
    segments_per_pole = 4,
    flux_scale = 1.0,
    end_effect_factor = 1.0,
    rotor_count = 2,
    field_model = :ideal_sheet,
    magnet_width_ratio = 1.0,
    field_eval_offset = 0.0,
)
    k_halbach = pi / tau_p
    segment_count = _smooth_max(segments_per_pole, 1.0)
    segmentation_factor = _halbach_segmentation_factor(segment_count, magnet_width_ratio, field_model)
    magnet_thickness_factor = 1 - exp(-k_halbach * h_m / mu_r)
    gap_factor = exp(-k_halbach * (len_ag + field_eval_offset))

    return rotor_count * flux_scale * end_effect_factor * B_r * segmentation_factor * magnet_thickness_factor * gap_factor
end

function PMSG_axial_Halbach(
    r_in,
    r_out,
    h_s,
    tau_p,
    h_m,
    h_ys,
    h_yr,
    machine_rating,
    shaft_rpm,
    Torque,
    b_st,
    d_s,
    t_ws,
    n_r,
    n_s,
    b_r,
    d_r,
    t_wr,
    D_shaft,
    rho_Fe,
    rho_Copper,
    rho_Fes,
    rho_PM;
    len_ag = 0.00075 * (r_in + r_out),
    B_r = 1.2,
    halbach_flux_boost = 1.0,          # optional calibration multiplier, unity by default
    halbach_field_model = :ideal_sheet,# :ideal_sheet, :finite_width_harmonic, or geometry-based :segmented_cuboid
    halbach_segments_per_pole = 4,     # magnetization steps per pole; larger approaches continuous Halbach
    halbach_magnets_total = nothing,   # explicit block count for :segmented_cuboid; defaults to 2*p*segments_per_pole
    halbach_magnetization_rotation = -1.0, # standard +z,-t,-z,+t sequence when negative
    halbach_coil_quadrature_order = 6,
    halbach_rotor_samples = 48,
    halbach_field_eval_offset = 0.0,   # extra distance from mechanical air gap to the winding/field evaluation plane [m]
    halbach_end_effect_factor = 1.0,   # finite-radius/end-effect derating when known from FEM or tests
    halbach_weak_side_fraction = 0.05, # residual weak-side flux crossing rotor back iron
    backiron_fraction = 0.5,           # fraction of rotor back-iron thickness retained
    turns_per_phase = nothing,         # explicit stator turns per phase for hand-wound/fractional-slot machines
    effective_flux_area = nothing,     # explicit flux area per pole/coil side [m^2]
    winding_factor = nothing,          # explicit winding factor; defaults to the distributed-winding estimate
    phase_resistance = nothing,        # measured or externally calculated phase resistance [Ohm]
    phase_inductance = nothing,        # measured or externally calculated phase inductance [H]
    wire_gauge_awg = nothing,          # bare round-wire AWG; used only when conductor_area is omitted
    conductor_area = nothing,          # bare area of one conductor/path [m^2]; overrides wire_gauge_awg
    wire_outer_diameter = nothing,     # insulated winding-wire diameter [m]; required for physical turn packing
    turns_per_layer = nothing,         # in-plane nested turns before starting a new axial layer
    turn_pitch = nothing,              # in-plane center spacing [m]; defaults to wire_outer_diameter
    layer_pitch = nothing,             # axial layer center spacing [m]; defaults to wire_outer_diameter
    winding_geometry_reference = :mean_turn, # :mean_turn or :inner_support
    winding_support_clearance = 0.0,   # bobbin/insulation distance to the first wire surface [m]
    mean_turn_length = nothing,         # full conductor length of one turn [m]
    coil_inner_radius = nothing,       # optional trapezoid inner radius for mean-turn calculation [m]
    coil_outer_radius = nothing,       # optional trapezoid outer radius for mean-turn calculation [m]
    coil_span_angle = nothing,         # optional trapezoid angular span for mean-turn calculation [rad]
    phase_lead_length = 0.0,           # additional series conductor length per phase path [m]
    phase_joint_resistance = 0.0,      # terminals/splices in one phase path [Ohm]
    parallel_paths = 1.0,              # equal electrical paths in parallel within one phase
    copper_resistivity_20c = nothing,  # when supplied, temperature-corrected instead of using resist_Cu
    copper_temperature_c = 20.0,
    copper_temperature_coefficient = 0.00393,
    inductance_model = :legacy_distributed, # :concentrated_coreless or physical :coreless_filament
    turns_per_coil = nothing,          # required for :concentrated_coreless
    coils_in_series_per_phase = nothing,# required for :concentrated_coreless
    inductance_coil_area = nothing,    # linked area of one concentrated coil [m^2]
    coil_mutual_coupling = 0.0,        # equal-pair mutual/self sensitivity; zero means uncoupled coils
    phase_leakage_inductance = 0.0,    # explicit phase leakage addition for concentrated model [H]
    phase_lead_inductance = 0.0,       # explicit lead contribution for :coreless_filament [H]
    inductance_path_subdivisions = 12,
    phase_coil_offset = 0.0,
    flux_linkage_factor = 1.0,         # spatial/3-D linked-flux factor; scales phi/E but not reported B_g
    airgap_flux_density = nothing,     # explicit fundamental winding-plane B [T]; overrides the Halbach field estimate
    phase_flux_linkage = nothing,      # explicit peak phase linkage [Wb-turn]; highest-priority flux input
    E = 2.0e11,
    P_Fe0e = 1.0,
    P_Fe0h = 4.0,
    alpha_p = pi / 2 * 0.7,
    b_s_tau_s = 0.45,
    b_so = 0.004,
    cofi = 0.85,
    h_i = 0.001,
    h_w = 0.005,
    k_fes = 0.9,
    k_fills = 0.65,
    m = 3.0,
    mu_0 = pi * 4e-7,
    mu_r = 1.06,
    phi = 90 * 2 * pi / 360.0,
    q1 = 1.0,
    ratio_mw2pp = 0.7,
    resist_Cu = 1.8 * 10^(-8) * 1.4,
    sigma = 40.0e3,
    gravity = 9.81,
    y_tau_p = 1.0,
    main_shaft_cm = [0.0, 0.0, 0.0],
    main_shaft_length = 2.0,
    v_poisson = 0.3,
    continuous = false,
    convergefaster = false,
    dual_rotor = true,
)
    R_sh = 0.5 * D_shaft
    Rm = 0.5 * (r_in + r_out)
    dr = r_out - r_in
    dr_eff = _smooth_max(dr, 1.0e-9)
    area_ag = pi * (r_out^2 - r_in^2)

    K_rad = dr / (r_out + r_in)

    if continuous
        p = pi * (r_in + r_out) / (2 * tau_p)
    else
        p = round(pi * (r_in + r_out) / (2 * tau_p))
    end

    f = shaft_rpm * p / 60.0
    S = 2 * p * q1 * m
    N_conductors = S * 2
    N_s_auto = N_conductors / (2 * m)
    N_s = turns_per_phase === nothing ? N_s_auto : turns_per_phase
    tau_s = 2 * pi * Rm / S
    b_s = b_s_tau_s * tau_s
    b_t = tau_s - b_s
    Slot_aspect_ratio = h_s / b_s

    ahm = len_ag + h_m / mu_r
    ba = b_so / (2 * ahm)
    gamma = 4 / pi * (ba * atan(ba) - log(sqrt(1 + ba^2)))
    k_C = tau_s / (tau_s - gamma * ahm)
    g_eff = k_C * ahm

    om_m = 2 * pi * shaft_rpm / 60.0
    om_e = p * om_m

    h_yr_eff = h_yr * backiron_fraction

    h_yr_safe = _smooth_max(h_yr_eff, 1.0e-6)
    rotor_count = dual_rotor ? 2 : 1
    spatial_halbach = if halbach_field_model == :segmented_cuboid
        if dual_rotor
            throw(ArgumentError(":segmented_cuboid currently represents a single free-space rotor; dual_rotor must be false"))
        end
        if coil_inner_radius === nothing || coil_outer_radius === nothing || coil_span_angle === nothing ||
           turns_per_coil === nothing || coils_in_series_per_phase === nothing
            throw(ArgumentError(":segmented_cuboid requires coil geometry, turns_per_coil, and coils_in_series_per_phase"))
        end
        magnets_total_use = halbach_magnets_total === nothing ? 2 * p * halbach_segments_per_pole : halbach_magnets_total
        magnet_width = ratio_mw2pp * tau_p / halbach_segments_per_pole
        segmented_halbach_winding_properties(;
            pole_pairs = p,
            magnet_inner_radius = r_in,
            magnet_outer_radius = r_out,
            magnet_thickness = h_m,
            magnet_tangential_width = magnet_width,
            magnets_total = magnets_total_use,
            B_r = B_r * halbach_flux_boost * halbach_end_effect_factor,
            winding_plane_offset = len_ag + halbach_field_eval_offset,
            coil_inner_radius,
            coil_outer_radius,
            coil_span_angle,
            turns_per_coil,
            coils_in_series_per_phase,
            phase_coil_offset,
            wire_outer_diameter,
            turns_per_layer,
            turn_pitch,
            layer_pitch,
            winding_geometry_reference,
            winding_support_clearance,
            phase_count = m,
            magnetization_rotation = halbach_magnetization_rotation,
            coil_quadrature_order = halbach_coil_quadrature_order,
            rotor_samples = halbach_rotor_samples,
        )
    else
        nothing
    end
    B_pm1 = spatial_halbach === nothing ? halbach_fundamental_flux_density(
        B_r,
        h_m,
        len_ag,
        tau_p;
        mu_r,
        segments_per_pole = halbach_segments_per_pole,
        flux_scale = halbach_flux_boost,
        end_effect_factor = halbach_end_effect_factor,
        rotor_count,
        field_model = halbach_field_model,
        magnet_width_ratio = ratio_mw2pp,
        field_eval_offset = halbach_field_eval_offset,
    ) : spatial_halbach.B_fundamental
    if airgap_flux_density !== nothing && airgap_flux_density <= 0
        throw(ArgumentError("airgap_flux_density must be positive"))
    end
    if phase_flux_linkage !== nothing && phase_flux_linkage <= 0
        throw(ArgumentError("phase_flux_linkage must be positive"))
    end
    B_g = airgap_flux_density === nothing ? B_pm1 : airgap_flux_density
    l_u = k_fes * dr
    l_e = dr
    b_m = ratio_mw2pp * tau_p
    B_symax = B_g * b_m * l_e / (2 * h_ys * l_u)
    B_rymax = halbach_weak_side_fraction * B_g * b_m * l_e / (2 * h_yr_safe * dr_eff)
    B_tmax = B_g * tau_s / b_t

    k_wd_auto = sin(pi / 6) / q1 / sin(pi / 6 / q1)
    k_wd = winding_factor === nothing ? k_wd_auto : winding_factor

    l_turn_legacy = 2 * dr + 2 * tau_p
    L_t = l_turn_legacy
    geometry_inputs = (coil_inner_radius, coil_outer_radius, coil_span_angle)
    geometry_input_count = count(x -> x !== nothing, geometry_inputs)
    if geometry_input_count != 0 && geometry_input_count != 3
        throw(ArgumentError("coil_inner_radius, coil_outer_radius, and coil_span_angle must be provided together"))
    end
    if geometry_input_count == 3 && (coil_outer_radius <= coil_inner_radius || coil_inner_radius < 0 || coil_span_angle <= 0)
        throw(ArgumentError("coil support geometry must have 0 <= inner radius < outer radius and positive span angle"))
    end
    if mean_turn_length !== nothing && mean_turn_length <= 0
        throw(ArgumentError("mean_turn_length must be positive"))
    end
    if phase_lead_length < 0 || phase_joint_resistance < 0 || parallel_paths <= 0
        throw(ArgumentError("lead length/joint resistance must be nonnegative and parallel_paths must be positive"))
    end
    if conductor_area !== nothing && conductor_area <= 0
        throw(ArgumentError("conductor_area must be positive"))
    end
    if wire_gauge_awg !== nothing && wire_gauge_awg < 0
        throw(ArgumentError("wire_gauge_awg must be nonnegative"))
    end
    if wire_outer_diameter !== nothing && wire_outer_diameter <= 0
        throw(ArgumentError("wire_outer_diameter must be positive"))
    end
    if winding_support_clearance < 0
        throw(ArgumentError("winding_support_clearance must be nonnegative"))
    end
    if copper_resistivity_20c !== nothing && copper_resistivity_20c <= 0
        throw(ArgumentError("copper_resistivity_20c must be positive"))
    end

    l_turn_physical = if mean_turn_length !== nothing
        mean_turn_length
    elseif geometry_input_count == 3
        2 * (coil_outer_radius - coil_inner_radius) + coil_span_angle * (coil_inner_radius + coil_outer_radius)
    else
        l_turn_legacy
    end
    supplied_conductor_area = conductor_area !== nothing ? conductor_area :
                              (wire_gauge_awg !== nothing ? _awg_conductor_area(wire_gauge_awg) : nothing)
    physical_winding_path = supplied_conductor_area !== nothing || mean_turn_length !== nothing || geometry_input_count == 3

    A_s = b_s * (h_s - h_w) * q1 * p
    A_scalc = b_s * 1000 * (h_s - h_w) * 1000 * q1 * p
    A_Cus_legacy = A_s * k_fills / N_s
    A_Cuscalc_legacy = A_scalc * k_fills / N_s
    packed_turns = if turns_per_layer !== nothing && turns_per_coil !== nothing && coils_in_series_per_phase !== nothing &&
                      geometry_input_count == 3 && wire_outer_diameter !== nothing
        _wedge_turn_geometries(
            coil_inner_radius,
            coil_outer_radius,
            coil_span_angle,
            turns_per_coil,
            turns_per_layer,
            turn_pitch === nothing ? wire_outer_diameter : turn_pitch,
            layer_pitch === nothing ? wire_outer_diameter : layer_pitch,
            winding_geometry_reference,
            winding_support_clearance,
        )
    else
        nothing
    end
    if packed_turns !== nothing
        _validate_interleaved_winding_clearance(
            packed_turns, m, coils_in_series_per_phase, wire_outer_diameter,
        )
    end
    phase_path_length = if packed_turns === nothing
        N_s * l_turn_physical + phase_lead_length
    else
        coils_in_series_per_phase * sum(_wedge_turn_length(turn) for turn in packed_turns) + phase_lead_length
    end
    l_Cus = physical_winding_path ? phase_path_length * parallel_paths : 2 * N_s * l_turn_legacy
    A_Cus = supplied_conductor_area === nothing ? A_Cus_legacy : supplied_conductor_area
    A_Cuscalc = supplied_conductor_area === nothing ? A_Cuscalc_legacy : supplied_conductor_area * 1.0e6
    winding_resistivity = copper_resistivity_20c === nothing ? resist_Cu :
        copper_resistivity_at_temperature(
            copper_resistivity_20c,
            copper_temperature_c;
            alpha_20 = copper_temperature_coefficient,
        )
    R_s_calc = physical_winding_path ?
        winding_resistivity * phase_path_length / (A_Cus * parallel_paths) + phase_joint_resistance :
        l_Cus * winding_resistivity / A_Cus
    R_s = phase_resistance === nothing ? R_s_calc : phase_resistance

    L_m_legacy = mu_0 * k_wd^2 * N_s^2 * area_ag / (g_eff * p)
    L_ssigmas = 2 * mu_0 * N_s^2 / p / q1 * dr * ((h_s - h_w) / (3 * b_s) + h_w / b_so)
    L_ssigmaew = 2 * mu_0 * N_s^2 / p / q1 * dr * 0.34 * len_ag * (l_e - 0.64 * tau_p * y_tau_p) / dr_eff
    L_ssigmag = 2 * mu_0 * N_s^2 / p / q1 * dr * (5 * (len_ag * k_C / b_so) / (5 + 4 * (len_ag * k_C / b_so)))
    L_ssigma = L_ssigmas + L_ssigmaew + L_ssigmag
    L_s_legacy = L_m_legacy + L_ssigma
    L_m = L_m_legacy
    L_s_calc = L_s_legacy
    if phase_inductance === nothing
        if inductance_model == :legacy_distributed
            L_s = L_s_legacy
        elseif inductance_model == :concentrated_coreless
            if turns_per_coil === nothing || coils_in_series_per_phase === nothing || inductance_coil_area === nothing
                throw(ArgumentError("turns_per_coil, coils_in_series_per_phase, and inductance_coil_area are required for :concentrated_coreless"))
            end
            if turns_per_coil <= 0 || coils_in_series_per_phase <= 0 || inductance_coil_area <= 0
                throw(ArgumentError("concentrated-coreless winding inputs must be positive"))
            end
            if abs(turns_per_coil * coils_in_series_per_phase - N_s) > 1.0e-8 * max(abs(N_s), 1.0)
                throw(ArgumentError("turns_per_coil * coils_in_series_per_phase must equal turns_per_phase"))
            end
            if coil_mutual_coupling < -inv(coils_in_series_per_phase - 1 + 1.0e-12)
                throw(ArgumentError("coil_mutual_coupling makes the phase magnetizing inductance negative"))
            end
            if phase_leakage_inductance < 0
                throw(ArgumentError("phase_leakage_inductance must be nonnegative"))
            end
            coil_self_inductance = mu_0 * k_wd^2 * turns_per_coil^2 * inductance_coil_area / g_eff
            mutual_factor = 1 + coil_mutual_coupling * (coils_in_series_per_phase - 1)
            L_m = coils_in_series_per_phase * coil_self_inductance * mutual_factor
            L_s_calc = L_m + phase_leakage_inductance
            L_s = L_s_calc
        elseif inductance_model == :coreless_filament
            if turns_per_coil === nothing || coils_in_series_per_phase === nothing || turns_per_layer === nothing ||
               wire_outer_diameter === nothing || geometry_input_count != 3
                throw(ArgumentError(":coreless_filament requires coil geometry, turns/layer layout, and wire_outer_diameter"))
            end
            if abs(turns_per_coil * coils_in_series_per_phase - N_s) > 1.0e-8 * max(abs(N_s), 1.0)
                throw(ArgumentError("turns_per_coil * coils_in_series_per_phase must equal turns_per_phase"))
            end
            inductance_result = coreless_winding_inductance(;
                coil_inner_radius,
                coil_outer_radius,
                coil_span_angle,
                turns_per_coil,
                coils_in_series_per_phase,
                wire_diameter = wire_outer_diameter,
                conductor_diameter = supplied_conductor_area === nothing ?
                    wire_outer_diameter : 2 * sqrt(supplied_conductor_area / pi),
                turns_per_layer,
                turn_pitch = turn_pitch === nothing ? wire_outer_diameter : turn_pitch,
                layer_pitch = layer_pitch === nothing ? wire_outer_diameter : layer_pitch,
                winding_geometry_reference,
                winding_support_clearance,
                path_subdivisions = inductance_path_subdivisions,
                phase_count = m,
                phase_lead_inductance,
                mu_0,
            )
            L_m = inductance_result.winding_inductance
            L_s_calc = inductance_result.phase_inductance
            L_s = L_s_calc
        else
            throw(ArgumentError("inductance_model must be :legacy_distributed, :concentrated_coreless, or :coreless_filament"))
        end
    else
        L_s = phase_inductance
    end

    flux_area_width_factor = halbach_field_model == :finite_width_harmonic ? 1.0 : ratio_mw2pp
    flux_area = effective_flux_area === nothing ? area_ag / (2 * p) * flux_area_width_factor : effective_flux_area
    lambda_phase = if phase_flux_linkage === nothing
        if flux_linkage_factor <= 0
            throw(ArgumentError("flux_linkage_factor must be positive"))
        end
        if spatial_halbach === nothing
            N_s * k_wd * B_g * flux_area * flux_linkage_factor
        else
            spatial_halbach.phase_flux_linkage_fundamental * flux_linkage_factor
        end
    else
        phase_flux_linkage
    end
    E_p = 4.44 * f * lambda_phase

    Z = machine_rating / (m * E_p)
    if convergefaster
        G = _smooth_abs((1.1 * E_p)^4 - (1 / 9) * (machine_rating * om_e * L_s)^2)
    else
        G = _smooth_max(E_p^2 - (om_e * L_s * Z)^2, 1.0e-6)
    end

    if convergefaster
        I_s = sqrt(2 * _smooth_abs((E_p * 1.1)^2 - sqrt(G)) / (om_e * L_s)^2)
    else
        I_s = sqrt(Z^2 + ((E_p - sqrt(G)) / (om_e * L_s))^2)
    end
    J_s = I_s / (A_Cuscalc * (physical_winding_path ? parallel_paths : 1.0))
    A_1 = 6 * N_s * I_s / (pi * 2 * Rm)

    B_smax = sqrt(2) * I_s * mu_0 / g_eff

    V_Cus = m * l_Cus * A_Cus
    V_Fest = dr * 2 * p * q1 * m * b_t * h_s
    V_Fesy = pi * (r_out^2 - r_in^2) * h_ys
    V_Fery = pi * (r_out^2 - r_in^2) * h_yr_eff
    Copper = V_Cus * rho_Copper
    M_Fest = V_Fest * rho_Fe
    M_Fesy = V_Fesy * rho_Fe
    M_Fery = V_Fery * rho_Fe
    Iron = M_Fest + M_Fesy + M_Fery

    mass_PM = area_ag * h_m * ratio_mw2pp * rho_PM * (dual_rotor ? 2 : 1)

    K_R = 1.2
    I_snom = machine_rating / (m * E_p * cofi)
    P_Cu = m * I_snom^2 * R_s * K_R

    P_Hyys = M_Fesy * (B_symax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60))
    P_Ftys = M_Fesy * (B_symax / 1.5)^2 * (P_Fe0e * (om_e / (2 * pi * 60))^2)
    P_Fesynom = P_Hyys + P_Ftys

    P_Hyd = M_Fest * (B_tmax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60))
    P_Ftd = M_Fest * (B_tmax / 1.5)^2 * (P_Fe0e * (om_e / (2 * pi * 60))^2)
    P_Festnom = P_Hyd + P_Ftd

    P_Hyyr = M_Fery * (B_rymax / 1.5)^2 * (P_Fe0h * om_e / (2 * pi * 60))
    P_Ftyr = M_Fery * (B_rymax / 1.5)^2 * (P_Fe0e * (om_e / (2 * pi * 60))^2)
    P_Ferynom = P_Hyyr + P_Ftyr

    P_ad = 0.2 * (P_Hyys + P_Ftys + P_Hyd + P_Ftd + P_Hyyr + P_Ftyr)
    pFtm = 300.0
    magnet_loss_area = 2 * p * b_m * dr * rotor_count
    P_Ftm = pFtm * magnet_loss_area * (f / 60)^2 * (B_g / _smooth_max(B_r, 1.0e-9))^2

    Losses = P_Cu + P_Festnom + P_Fesynom + P_ad + P_Ftm + P_Ferynom
    gen_eff = machine_rating / (machine_rating + Losses)

    q3 = B_g^2 / (2 * mu_0)

    r_plate_inner = _smooth_max(R_sh, r_in)
    u_ar = plate_deflection_uniform(q3, r_plate_inner, r_out + h_yr_eff, h_yr_safe, E; v = v_poisson)
    u_as = plate_deflection_uniform(q3, R_sh, r_in, h_ys, E; v = v_poisson)

    y_ar = u_ar
    y_as = u_as

    G_mod = E / (2 * (1 + v_poisson))
    J_rotor = 0.5 * pi * ((r_out + h_yr_eff)^4 - r_plate_inner^4)
    J_stator = 0.5 * pi * ((r_out + h_ys)^4 - (R_sh)^4)
    theta_r = Torque * h_yr_safe / (J_rotor * G_mod)
    theta_s = Torque * h_ys / (J_stator * G_mod)
    z_ar = theta_r * r_out
    z_as = theta_s * r_out

    u_allow_r = r_out / 10000
    u_allow_s = r_out / 10000
    y_allow = dr * 0.02
    z_allow_s = 0.05 * 2 * pi * r_out / 360
    z_allow_r = 0.05 * 2 * pi * r_out / 360
    b_allow_s = 2 * pi * R_sh / _smooth_max(1.0, n_s)
    b_allow_r = 2 * pi * R_sh / _smooth_max(1.0, n_r)

    a_r = (b_r * d_r) - ((b_r - 2 * t_wr) * (d_r - 2 * t_wr))
    a_s = (b_st * d_s) - ((b_st - 2 * t_ws) * (d_s - 2 * t_ws))
    mass_arms = (n_r * a_r + n_s * a_s) * dr * rho_Fes
    mass_rotor_disc = rho_Fes * pi * ((r_out + h_yr_eff)^2 - r_plate_inner^2) * h_yr_eff
    mass_stator_disc = rho_Fes * pi * ((r_out + h_ys)^2 - (R_sh)^2) * h_ys
    Structural_mass = mass_arms + mass_rotor_disc + mass_stator_disc

    Stator = M_Fesy + M_Fest + Copper
    Rotor = M_Fery + mass_PM
    Mass = Stator + Rotor + Structural_mass

    TC1 = Torque / (2 * pi * sigma)
    TC2r = Rm^2 * dr
    TC2s = Rm^2 * dr

    R_out = r_out + h_m + _smooth_max(h_yr_eff, h_ys)
    len_axial_total = h_yr_eff + h_ys + (dual_rotor ? 2 * h_m : h_m) + 2 * len_ag

    I = zeros(Real, 3)
    I[1] = 0.5 * Mass * R_out^2
    I[2] = 0.25 * Mass * R_out^2 + (1 / 12) * Mass * len_axial_total^2
    I[3] = I[2]
    cm = zeros(Real, 3)
    cm[1] = main_shaft_cm[1] + main_shaft_length / 2.0 + len_axial_total / 2
    cm[2] = main_shaft_cm[2]
    cm[3] = main_shaft_cm[3]

    R_1 = r_in

    return B_symax,
    B_tmax,
    B_rymax,
    B_smax,
    B_pm1,
    B_g,
    N_s,
    b_s,
    b_t,
    A_Cuscalc,
    b_m,
    p,
    E_p,
    f,
    I_s,
    R_s,
    L_s,
    A_1,
    J_s,
    Losses,
    K_rad,
    gen_eff,
    S,
    Slot_aspect_ratio,
    Copper,
    Iron,
    u_ar,
    y_ar,
    z_ar,
    u_as,
    y_as,
    z_as,
    u_allow_r,
    u_allow_s,
    y_allow,
    z_allow_s,
    z_allow_r,
    b_allow_s,
    b_allow_r,
    TC1,
    TC2r,
    TC2s,
    R_out,
    Structural_mass,
    Mass,
    mass_PM,
    cm,
    I,
    R_1
end
