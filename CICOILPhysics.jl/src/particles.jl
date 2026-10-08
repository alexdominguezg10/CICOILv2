# OilParticles — Structure-of-Arrays (SoA) layout for particle state.
#
# Uses AbstractVector/AbstractMatrix so the same struct works for both
# CPU (Vector/Matrix) and GPU (CuVector/CuMatrix) backends.

struct OilParticles{V <: AbstractVector{Float32}, M <: AbstractMatrix{Float32},
                    VU <: AbstractVector{UInt8}}
    # ── position / geometry ──────────────────────────────────────────────────
    z               ::V    # depth [m], ≤0 (0 = surface)
    diameter        ::V    # droplet diameter [m], 0 for surface spillets
    age_seconds     ::V    # particle age [s]
    is_dissolved    ::VU   # 1 = dissolved TAMOC far-field particle

    # ── per-particle mixture scalars ─────────────────────────────────────────
    oil_molar_mass    ::V   # mixture Mw  [kg/mol]
    oil_density       ::V   # mixture rho [kg/m³]
    spillet_thickness ::V   # film thickness [m]
    water_fraction    ::V   # emulsion water content [0–1]
    kinematic_viscosity::V  # [cSt]
    max_water         ::V   # per-particle emulsification cap [0–1]

    # ── mass sinks — scalar totals ───────────────────────────────────────────
    mass_oil                   ::V   # remaining oil mass [kg]
    mass_evaporated            ::V
    mass_dispersed             ::V
    mass_biodegraded           ::V   # total (oil + dissolved)
    mass_biodegraded_from_oil  ::V
    mass_biodegraded_from_water::V
    mass_photooxidized         ::V
    mass_op_dissolved          ::V   # dissolved OP pool [kg]
    mass_op_degraded           ::V

    # ── fractions (diagnostics) ──────────────────────────────────────────────
    fraction_evaporated        ::V
    fraction_biodegraded_from_oil::V
    fraction_photooxidized     ::V
    fraction_op_dissolved      ::V
    fraction_op_degraded       ::V

    # ── 17-component mass balance (N × 17) ───────────────────────────────────
    mass_components    ::M   # current component masses [kg]
    mass_evap_comp     ::M   # accumulated evaporated per component
    mass_biodeg_comp   ::M   # accumulated biodegraded (oil phase)
    mass_photoox_comp  ::M   # accumulated photooxidized
    mass_op_diss_comp  ::M   # dissolved OP pool per component
    mass_op_deg_comp   ::M   # degraded OP per component
    mass_dissolved_comp::M   # dissolved-phase mass (TAMOC particles)
end

"""
    OilParticles(N; initial_mass_kg=185.6f0, oil_props=nothing, b=backend())

Allocate an `OilParticles` for `N` particles on the active backend.
"""
function OilParticles(N::Int;
    initial_mass_kg::Float32 = 185.6f0,
    oil_props::Union{OilProps, Nothing} = nothing,
    b::AbstractBackend = backend(),
)
    init_comp_host = if isnothing(oil_props)
        fill(initial_mass_kg / 17.0f0, N, 17)
    else
        Float32.(ones(N, 1) * collect(Float32, oil_props.oil_mass_frac .* initial_mass_kg)')
    end

    az(N) = alloc_vector(b, Float32, N)
    am(N, M) = alloc_matrix(b, Float32, N, M)
    fv(val, N) = fill_vector(b, val, N)

    OilParticles(
        az(N),                # z
        az(N),                # diameter
        az(N),                # age_seconds
        alloc_vector(b, UInt8, N),  # is_dissolved

        fv(200.0f0, N),       # oil_molar_mass
        fv(880.0f0, N),       # oil_density
        fv(0.01f0, N),        # spillet_thickness
        az(N),                # water_fraction
        fv(50.0f0, N),        # kinematic_viscosity
        fv(EMUL_MAX_WATER_DEFAULT, N),  # max_water

        fv(initial_mass_kg, N),  # mass_oil
        az(N),                # mass_evaporated
        az(N),                # mass_dispersed
        az(N),                # mass_biodegraded
        az(N),                # mass_biodegraded_from_oil
        az(N),                # mass_biodegraded_from_water
        az(N),                # mass_photooxidized
        az(N),                # mass_op_dissolved
        az(N),                # mass_op_degraded

        az(N),                # fraction_evaporated
        az(N),                # fraction_biodegraded_from_oil
        az(N),                # fraction_photooxidized
        az(N),                # fraction_op_dissolved
        az(N),                # fraction_op_degraded

        to_backend(b, init_comp_host),  # mass_components (N, 17)
        am(N, 17),            # mass_evap_comp
        am(N, 17),            # mass_biodeg_comp
        am(N, 17),            # mass_photoox_comp
        am(N, 17),            # mass_op_diss_comp
        am(N, 17),            # mass_op_deg_comp
        am(N, 17),            # mass_dissolved_comp
    )
end

struct Environment{V <: AbstractVector{Float32}}
    sea_water_temperature ::V   # K
    x_sea_water_velocity  ::V   # m/s
    y_sea_water_velocity  ::V
    x_wind                ::V   # m/s
    y_wind                ::V
    sea_surface_wave_significant_height::V  # Hs [m]
    sea_surface_wave_period_at_variance_spectral_density_maximum::V  # Tp [s]
    uv_irradiance         ::Float32  # W/m², scalar
end

"""
    Environment(env_dict; uv_irradiance=0.0f0, b=backend())

Build Environment from a Dict{String,AbstractArray}.
"""
function Environment(env::Dict{String, <:AbstractArray};
    uv_irradiance::Float32=0.0f0,
    b::AbstractBackend = backend(),
)
    tb(x) = to_backend(b, Float32.(x))
    Environment(
        tb(env["sea_water_temperature"]),
        tb(env["x_sea_water_velocity"]),
        tb(env["y_sea_water_velocity"]),
        tb(env["x_wind"]),
        tb(env["y_wind"]),
        tb(env["sea_surface_wave_significant_height"]),
        tb(env["sea_surface_wave_period_at_variance_spectral_density_maximum"]),
        Float32(uv_irradiance),
    )
end

struct WeatheringConfig
    evaporation    ::Bool
    emulsification ::Bool
    biodegradation ::Bool
    photooxidation ::Bool
end

WeatheringConfig(; evaporation=true, emulsification=true,
                   biodegradation=false, photooxidation=false) =
    WeatheringConfig(evaporation, emulsification, biodegradation, photooxidation)

nparticles(p::OilParticles) = length(p.z)
