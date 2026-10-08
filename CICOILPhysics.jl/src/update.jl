# update_weathering! — master dispatch for one simulation timestep.

"""
    update_weathering!(p, env, props, config, transfer_coeff, dt; b=backend())

Apply one timestep of weathering physics. Backend `b` selects CPU or GPU.
"""
function update_weathering!(
    p             ::OilParticles,
    env           ::Environment,
    props         ::OilProps,
    config        ::WeatheringConfig,
    transfer_coeff::AbstractVector{Float32},
    dt            ::Float32;
    b             ::AbstractBackend = backend(),
)
    config.evaporation    && launch_evaporation!(p, env, props, transfer_coeff, dt; b)
    config.emulsification && launch_emulsification!(p, env, dt; b)
    config.biodegradation && launch_biodegradation!(p, env, dt; b)
    config.photooxidation && launch_photooxidation!(p, env, dt; b)
    sync_backend!(b)
    return
end

# ── Transfer coefficient (always computed on CPU) ────────────────────────────

function compute_transfer_coeff(
    x_wind ::AbstractVector{Float32},
    y_wind ::AbstractVector{Float32},
    oil_molar_mass::AbstractVector{Float32},
    sea_water_temp::AbstractVector{Float32},
    atm_pressure  ::Float32 = P_ATM;
    b             ::AbstractBackend = backend(),
)
    # Always compute on host (CPU)
    xw = to_host(b, x_wind)
    yw = to_host(b, y_wind)
    om = to_host(b, oil_molar_mass)
    st = to_host(b, sea_water_temp)

    N = length(xw)
    Kt = Vector{Float32}(undef, N)
    ν_air = 1.5f-5

    for i in 1:N
        Mw   = om[i]
        U10  = sqrt(xw[i]^2 + yw[i]^2)
        D_air = 4.0f-6 / sqrt(Mw / 30.0f0)
        Sc = ν_air / D_air
        Kt[i] = 0.0292f0 * max(U10, 0.5f0)^0.78f0 * (Sc / 0.6f0)^(-0.67f0)
    end

    return to_backend(b, Kt)
end
