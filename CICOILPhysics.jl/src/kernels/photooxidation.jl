# Photooxidation kernel — Beer-Lambert UV + aromatic group rate constants
#
# 1. OP removal (every timestep): dm_op = m_op × (1 - exp(-K_OP × q10 × dt_days))
# 2. Photooxidation (daytime, within UV depth): aromatic groups only (c=9..16)

function _photooxidation_body!(
    i::Int,
    mass_components  ::AbstractMatrix{Float32},
    mass_photoox_comp::AbstractMatrix{Float32},
    mass_op_diss_comp::AbstractMatrix{Float32},
    mass_op_deg_comp ::AbstractMatrix{Float32},
    mass_oil         ::AbstractVector{Float32},
    mass_photooxidized::AbstractVector{Float32},
    mass_op_dissolved::AbstractVector{Float32},
    mass_op_degraded ::AbstractVector{Float32},
    frac_photoox     ::AbstractVector{Float32},
    frac_op_diss     ::AbstractVector{Float32},
    frac_op_deg      ::AbstractVector{Float32},
    mass_evaporated  ::AbstractVector{Float32},
    mass_dispersed   ::AbstractVector{Float32},
    mass_biodeg      ::AbstractVector{Float32},
    z                ::AbstractVector{Float32},
    sea_water_temp   ::AbstractVector{Float32},
    k_photo          ::NTuple{8, Float32},
    I_UV0            ::Float32,
    dt               ::Float32,
)
    T_C = sea_water_temp[i] - 273.15f0
    q10 = BIODEG_Q10 ^ ((T_C - BIODEG_T_REF) / 10.0f0)
    dt_days = dt / 86400.0f0

    # ── 1. OP pool removal (unconditional) ────────────────────────────────────
    dm_op_deg_total = 0.0f0
    for c in 1:17
        m_op = mass_op_diss_comp[i, c]
        m_op <= 0.0f0 && continue

        k_op = K_OP_REMOVAL * q10
        dm   = m_op * (1.0f0 - exp(-k_op * dt_days))
        dm   = min(dm, m_op)

        mass_op_diss_comp[i, c] -= dm
        mass_op_deg_comp[i, c]  += dm
        dm_op_deg_total         += dm
    end
    mass_op_dissolved[i] -= dm_op_deg_total
    mass_op_degraded[i]  += dm_op_deg_total

    # ── 2. Photooxidation (daytime + within UV depth cutoff) ─────────────────
    if I_UV0 > 0.0f0 && z[i] >= -UV_DEPTH_CUTOFF
        depth   = abs(z[i])
        atten   = exp(-K_D_UV * depth)
        E_UV    = I_UV0 * atten * dt

        dm_phox_total = 0.0f0

        for c in 9:16
            m_c = mass_components[i, c]
            m_c <= 0.0f0 && continue

            k_c  = k_photo[c - 8]
            dm   = k_c * E_UV * m_c
            dm   = min(dm, m_c)
            dm   = max(dm, 0.0f0)

            dm_op  = dm * OP_DISSOLVED_FRACTION
            dm_res = dm - dm_op

            mass_components[i, c]    -= dm
            mass_photoox_comp[i, c]  += dm_res
            mass_op_diss_comp[i, c]  += dm_op

            dm_phox_total += dm
        end

        dm_phox_op  = dm_phox_total * OP_DISSOLVED_FRACTION
        dm_phox_res = dm_phox_total - dm_phox_op

        mass_oil[i]          -= dm_phox_total
        mass_photooxidized[i]+= dm_phox_res
        mass_op_dissolved[i] += dm_phox_op
    end

    # ── 3. Update diagnostic fractions ───────────────────────────────────────
    initial_mass = (mass_oil[i] + mass_evaporated[i] + mass_dispersed[i] +
                    mass_photooxidized[i] + mass_biodeg[i] +
                    mass_op_dissolved[i] + mass_op_degraded[i])

    if initial_mass > 0.0f0
        frac_photoox[i]  = mass_photooxidized[i]  / initial_mass
        frac_op_diss[i]  = mass_op_dissolved[i]   / initial_mass
        frac_op_deg[i]   = mass_op_degraded[i]    / initial_mass
    end
    return
end

function launch_photooxidation!(
    p     ::OilParticles,
    env   ::Environment,
    dt    ::Float32;
    b     ::AbstractBackend = backend(),
)
    N = nparticles(p)
    run_kernel!(b, _photooxidation_body!, N,
        p.mass_components, p.mass_photoox_comp,
        p.mass_op_diss_comp, p.mass_op_deg_comp,
        p.mass_oil, p.mass_photooxidized,
        p.mass_op_dissolved, p.mass_op_degraded,
        p.fraction_photooxidized, p.fraction_op_dissolved,
        p.fraction_op_degraded,
        p.mass_evaporated, p.mass_dispersed, p.mass_biodegraded,
        p.z, env.sea_water_temperature,
        K_PHOTO_AROM, env.uv_irradiance,
        dt,
    )
    return
end
