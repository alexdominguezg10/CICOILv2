# Evaporation kernel — Cox-Pv model (Cox 1923; Mackay & Matsugu 1973)
#
# For each surface particle (z==0, age < MAX_AGE_EVAP_DAYS):
#   Pv_c = P_atm * exp((1 - Tb_c/T) * exp(A0_c + A1_c*T + A2_c*T²))
#   x_c  = (m_c / Mw_c) / Σ(m_j / Mw_j)
#   ψ_c  = Kt/h * Mw_mix/ρ_mix * Pv_c * x_c / (R*T) * dt
#   dm_c = m_c * clamp(ψ_c, 0, 1)

function _evaporation_body!(
    i::Int,
    mass_components  ::AbstractMatrix{Float32},
    mass_evap_comp   ::AbstractMatrix{Float32},
    mass_oil         ::AbstractVector{Float32},
    mass_evaporated  ::AbstractVector{Float32},
    frac_evaporated  ::AbstractVector{Float32},
    mass_dispersed   ::AbstractVector{Float32},
    mass_biodeg      ::AbstractVector{Float32},
    mass_op_dissolved::AbstractVector{Float32},
    z                ::AbstractVector{Float32},
    age_seconds      ::AbstractVector{Float32},
    oil_molar_mass   ::AbstractVector{Float32},
    oil_density      ::AbstractVector{Float32},
    spillet_thickness::AbstractVector{Float32},
    transfer_coeff   ::AbstractVector{Float32},
    sea_water_temp   ::AbstractVector{Float32},
    comp_Mw ::NTuple{17, Float32},
    comp_Tb ::NTuple{17, Float32},
    Cox_A0  ::NTuple{17, Float32},
    Cox_A1  ::NTuple{17, Float32},
    Cox_A2  ::NTuple{17, Float32},
    dt      ::Float32,
)
    z[i] != 0.0f0                         && return
    age_seconds[i] > MAX_AGE_EVAP_DAYS * 86400.0f0 && return

    T  = sea_water_temp[i]
    Kt = transfer_coeff[i]
    h  = spillet_thickness[i]
    Mw = oil_molar_mass[i]
    ρ  = oil_density[i]

    total_moles = 0.0f0
    for c in 1:17
        total_moles += mass_components[i, c] / comp_Mw[c]
    end
    total_moles <= 0.0f0 && return

    m_evap_total = 0.0f0
    for c in 1:17
        m_c = mass_components[i, c]
        m_c <= 0.0f0 && continue

        x_c = (m_c / comp_Mw[c]) / total_moles
        Tb = comp_Tb[c]
        Pv = P_ATM * exp((1.0f0 - Tb / T) * exp(Cox_A0[c] + Cox_A1[c] * T + Cox_A2[c] * T * T))
        ψ = (Kt / h) * (Mw / (ρ * 1000.0f0)) * Pv * x_c / (R_GAS * T) * dt
        ψ = clamp(ψ, 0.0f0, 1.0f0)
        dm = m_c * ψ

        mass_components[i, c]   -= dm
        mass_evap_comp[i, c]    += dm
        m_evap_total            += dm
    end

    mass_oil[i]        -= m_evap_total
    mass_evaporated[i] += m_evap_total

    m_total = (mass_oil[i] + mass_evaporated[i] + mass_dispersed[i] +
               mass_biodeg[i] + mass_op_dissolved[i])
    if m_total > 0.0f0
        frac_evaporated[i] = mass_evaporated[i] / m_total
    end
    return
end

function launch_evaporation!(
    p     ::OilParticles,
    env   ::Environment,
    props ::OilProps,
    transfer_coeff::AbstractVector{Float32},
    dt    ::Float32;
    b     ::AbstractBackend = backend(),
)
    N = nparticles(p)
    run_kernel!(b, _evaporation_body!, N,
        p.mass_components, p.mass_evap_comp,
        p.mass_oil, p.mass_evaporated, p.fraction_evaporated,
        p.mass_dispersed, p.mass_biodegraded, p.mass_op_dissolved,
        p.z, p.age_seconds, p.oil_molar_mass, p.oil_density,
        p.spillet_thickness, transfer_coeff,
        env.sea_water_temperature,
        props.comp_Mw, props.comp_Tb,
        props.Cox_A0, props.Cox_A1, props.Cox_A2,
        dt,
    )
    return
end
