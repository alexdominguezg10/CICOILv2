# Python ↔ Julia bridge for juliacall
#
# Usage from Python:
#   from juliacall import Main as jl
#   jl.seval('using CICOILPhysics')
#   state = jl.CICOILPhysics.py_init_state(N_particles=3000, ...)
#   jl.CICOILPhysics.py_update_weathering(state, ...)

mutable struct SimState
    particles ::OilParticles
    props     ::OilProps
    N         ::Int
    b         ::AbstractBackend
end

function py_init_state(;
    N_particles    ::Int,
    oil_Mw         ::AbstractVector,
    oil_Tb         ::AbstractVector,
    oil_mass_frac  ::AbstractVector,
    initial_mass_kg::Real = 185.6,
)
    b = backend()
    props = OilProps(; Mw=oil_Mw, Tb=oil_Tb, oil_mass_frac=oil_mass_frac)
    particles = OilParticles(N_particles;
        initial_mass_kg = Float32(initial_mass_kg),
        oil_props = props,
        b = b,
    )
    SimState(particles, props, N_particles, b)
end

function py_update_weathering(
    state        ::SimState,
    z_py, diameter_py, age_seconds_py, is_dissolved_py,
    oil_molar_mass_py, oil_density_py, spillet_thickness_py,
    water_fraction_py, kinematic_visc_py, max_water_py,
    mass_oil_py, mass_evaporated_py, mass_dispersed_py,
    mass_biodeg_py, mass_biodeg_oil_py, mass_biodeg_water_py,
    mass_photoox_py, mass_op_diss_py, mass_op_deg_py,
    mass_components_py, mass_evap_comp_py, mass_biodeg_comp_py,
    mass_photoox_comp_py, mass_op_diss_comp_py, mass_op_deg_comp_py,
    mass_dissolved_comp_py,
    sea_water_temp_py, x_wind_py, y_wind_py, x_vel_py, y_vel_py,
    Hs_py, Tp_py,
    uv_irradiance  ::Real,
    do_evap ::Bool, do_emuls::Bool, do_biodeg::Bool, do_photoox::Bool,
    dt      ::Real,
)
    b = state.b
    p = state.particles
    tb(x) = to_backend(b, Float32.(x))

    # ── Upload particle scalars ──────────────────────────────────────────────
    copyto!(p.z,                  tb(z_py))
    copyto!(p.diameter,           tb(diameter_py))
    copyto!(p.age_seconds,        tb(age_seconds_py))
    copyto!(p.is_dissolved,       to_backend(b, UInt8.(is_dissolved_py)))
    copyto!(p.oil_molar_mass,     tb(oil_molar_mass_py))
    copyto!(p.oil_density,        tb(oil_density_py))
    copyto!(p.spillet_thickness,  tb(spillet_thickness_py))
    copyto!(p.water_fraction,     tb(water_fraction_py))
    copyto!(p.kinematic_viscosity,tb(kinematic_visc_py))
    copyto!(p.max_water,          tb(max_water_py))
    copyto!(p.mass_oil,           tb(mass_oil_py))
    copyto!(p.mass_evaporated,    tb(mass_evaporated_py))
    copyto!(p.mass_dispersed,     tb(mass_dispersed_py))
    copyto!(p.mass_biodegraded,   tb(mass_biodeg_py))
    copyto!(p.mass_biodegraded_from_oil,   tb(mass_biodeg_oil_py))
    copyto!(p.mass_biodegraded_from_water, tb(mass_biodeg_water_py))
    copyto!(p.mass_photooxidized, tb(mass_photoox_py))
    copyto!(p.mass_op_dissolved,  tb(mass_op_diss_py))
    copyto!(p.mass_op_degraded,   tb(mass_op_deg_py))

    copyto!(p.mass_components,    tb(mass_components_py))
    copyto!(p.mass_evap_comp,     tb(mass_evap_comp_py))
    copyto!(p.mass_biodeg_comp,   tb(mass_biodeg_comp_py))
    copyto!(p.mass_photoox_comp,  tb(mass_photoox_comp_py))
    copyto!(p.mass_op_diss_comp,  tb(mass_op_diss_comp_py))
    copyto!(p.mass_op_deg_comp,   tb(mass_op_deg_comp_py))
    copyto!(p.mass_dissolved_comp,tb(mass_dissolved_comp_py))

    # ── Build environment ────────────────────────────────────────────────────
    env = Environment(
        tb(sea_water_temp_py), tb(x_vel_py), tb(y_vel_py),
        tb(x_wind_py), tb(y_wind_py), tb(Hs_py), tb(Tp_py),
        Float32(uv_irradiance),
    )

    config = WeatheringConfig(
        evaporation=do_evap, emulsification=do_emuls,
        biodegradation=do_biodeg, photooxidation=do_photoox,
    )

    Kt = compute_transfer_coeff(
        env.x_wind, env.y_wind,
        p.oil_molar_mass, env.sea_water_temperature; b,
    )

    update_weathering!(p, env, state.props, config, Kt, Float32(dt); b)

    # ── Copy results back to Python/NumPy buffers ────────────────────────────
    h(x) = to_host(b, x)
    copyto!(mass_oil_py,          h(p.mass_oil))
    copyto!(mass_evaporated_py,   h(p.mass_evaporated))
    copyto!(mass_biodeg_py,       h(p.mass_biodegraded))
    copyto!(mass_biodeg_oil_py,   h(p.mass_biodegraded_from_oil))
    copyto!(mass_biodeg_water_py, h(p.mass_biodegraded_from_water))
    copyto!(mass_photoox_py,      h(p.mass_photooxidized))
    copyto!(mass_op_diss_py,      h(p.mass_op_dissolved))
    copyto!(mass_op_deg_py,       h(p.mass_op_degraded))
    copyto!(water_fraction_py,    h(p.water_fraction))
    copyto!(kinematic_visc_py,    h(p.kinematic_viscosity))
    copyto!(mass_components_py,   h(p.mass_components))
    copyto!(mass_evap_comp_py,    h(p.mass_evap_comp))
    copyto!(mass_biodeg_comp_py,  h(p.mass_biodeg_comp))
    copyto!(mass_photoox_comp_py, h(p.mass_photoox_comp))
    copyto!(mass_op_diss_comp_py, h(p.mass_op_diss_comp))
    copyto!(mass_op_deg_comp_py,  h(p.mass_op_deg_comp))

    return nothing
end
