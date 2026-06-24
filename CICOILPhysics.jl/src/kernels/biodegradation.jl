# Biodegradation kernel — per-component first-order decay (Adcroft et al. 2010)
#
# k_eff = K_BIODEG_ALL[c] × Q10^((T-T_REF)/10) × surf_factor × size_factor / 86400
# dm_c  = m_c × (1 - exp(-k_eff × dt))

function _biodegradation_body!(
    i::Int,
    mass_components      ::AbstractMatrix{Float32},
    mass_biodeg_comp     ::AbstractMatrix{Float32},
    mass_biodegraded_oil ::AbstractVector{Float32},
    mass_biodeg_total    ::AbstractVector{Float32},
    mass_dissolved_comp  ::AbstractMatrix{Float32},
    mass_biodeg_water    ::AbstractVector{Float32},
    mass_oil             ::AbstractVector{Float32},
    frac_biodeg_oil      ::AbstractVector{Float32},
    z                    ::AbstractVector{Float32},
    diameter             ::AbstractVector{Float32},
    is_dissolved         ::AbstractVector{UInt8},
    sea_water_temp       ::AbstractVector{Float32},
    k_all                ::NTuple{17, Float32},
    k_water              ::NTuple{17, Float32},
    dt                   ::Float32,
)
    T_C = sea_water_temp[i] - 273.15f0
    q10 = BIODEG_Q10 ^ ((T_C - BIODEG_T_REF) / 10.0f0)

    surf_factor = (z[i] == 0.0f0) ? BIODEG_SURFACE_FACTOR : 1.0f0

    d = diameter[i]
    d_eff = (d == 0.0f0) ? BIODEG_DROPLET_D_REF :
                            clamp(d, BIODEG_DROPLET_D_MIN, BIODEG_DROPLET_D_MAX)
    size_factor = (BIODEG_DROPLET_D_REF / d_eff) ^ BIODEG_DROPLET_EXPONENT

    dt_days = dt / 86400.0f0
    dm_oil_total = 0.0f0

    for c in 1:17
        m_c = mass_components[i, c]
        m_c <= 0.0f0 && continue

        k_eff = k_all[c] * q10 * surf_factor * size_factor
        dm    = m_c * (1.0f0 - exp(-k_eff * dt_days))
        dm    = min(dm, m_c)

        mass_components[i, c]   -= dm
        mass_biodeg_comp[i, c]  += dm
        dm_oil_total             += dm
    end

    mass_oil[i]              -= dm_oil_total
    mass_biodegraded_oil[i]  += dm_oil_total
    mass_biodeg_total[i]     += dm_oil_total

    m_remaining = mass_oil[i]
    m_total = m_remaining + mass_biodegraded_oil[i]
    if m_total > 0.0f0
        frac_biodeg_oil[i] = mass_biodegraded_oil[i] / m_total
    end

    is_dissolved[i] != 1 && return

    dm_water_total = 0.0f0
    for c in 1:17
        m_d = mass_dissolved_comp[i, c]
        m_d <= 0.0f0 && continue

        k_eff = k_water[c] * q10
        dm    = m_d * (1.0f0 - exp(-k_eff * dt_days))
        dm    = min(dm, m_d)

        mass_dissolved_comp[i, c] -= dm
        dm_water_total            += dm
    end

    mass_biodeg_water[i] += dm_water_total
    mass_biodeg_total[i] += dm_water_total
    return
end

function launch_biodegradation!(p::OilParticles, env::Environment, dt::Float32;
    b::AbstractBackend = backend(),
)
    N = nparticles(p)
    run_kernel!(b, _biodegradation_body!, N,
        p.mass_components, p.mass_biodeg_comp,
        p.mass_biodegraded_from_oil, p.mass_biodegraded,
        p.mass_dissolved_comp, p.mass_biodegraded_from_water,
        p.mass_oil, p.fraction_biodegraded_from_oil,
        p.z, p.diameter, p.is_dissolved,
        env.sea_water_temperature,
        K_BIODEG_ALL, K_BIODEG_WATER,
        dt,
    )
    return
end
