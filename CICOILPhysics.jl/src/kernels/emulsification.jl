# Emulsification kernel — Fingas & Fieldhouse (2004) viscosity-stability model
#
# For surface particles (z==0) with remaining oil:
#   k_emul_eff = k_emul × clip(ν_dry / NU_REF, MIN, MAX) × EMUL_RATE_SCALE
#   dW/dt      = k_emul_eff × (W_max - W)² × W_max
#   ν_emul     = ν_dry × exp(C_visc × W)

const MAYA_A_VISC = 0.00235f0
const MAYA_B_VISC = 3440.0f0
const C_MOONEY    = 2.5f0

function _emulsification_body!(
    i::Int,
    water_fraction   ::AbstractVector{Float32},
    kinematic_visc   ::AbstractVector{Float32},
    mass_oil         ::AbstractVector{Float32},
    z                ::AbstractVector{Float32},
    max_water        ::AbstractVector{Float32},
    sea_water_temp   ::AbstractVector{Float32},
    x_wind           ::AbstractVector{Float32},
    y_wind           ::AbstractVector{Float32},
    dt               ::Float32,
)
    z[i] != 0.0f0      && return
    mass_oil[i] <= 0.0f0 && return

    T  = sea_water_temp[i]
    Wm = max_water[i]
    W  = water_fraction[i]

    ν_dry = MAYA_A_VISC * exp(MAYA_B_VISC / T)
    Uwind = sqrt(x_wind[i]^2 + y_wind[i]^2)
    k_emul = 2.0f-6 * Uwind^2

    stability = clamp(ν_dry / EMUL_STABILITY_NU_REF, EMUL_STABILITY_MIN, EMUL_STABILITY_MAX)
    k_emul_eff = k_emul * stability * EMUL_RATE_SCALE

    dW = k_emul_eff * (Wm - W)^2 * Wm * dt
    W_new = min(W + dW, Wm)
    water_fraction[i] = W_new
    kinematic_visc[i] = ν_dry * exp(C_MOONEY * W_new)
    return
end

function launch_emulsification!(p::OilParticles, env::Environment, dt::Float32;
    b::AbstractBackend = backend(),
)
    N = nparticles(p)
    run_kernel!(b, _emulsification_body!, N,
        p.water_fraction, p.kinematic_viscosity,
        p.mass_oil, p.z, p.max_water,
        env.sea_water_temperature, env.x_wind, env.y_wind,
        dt,
    )
    return
end
