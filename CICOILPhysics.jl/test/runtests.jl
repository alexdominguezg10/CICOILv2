using Test

# Load CICOILPhysics from the parent directory
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using CICOILPhysics
using CICOILPhysics: K_BIODEG_ALL, K_PHOTO_AROM, BIODEG_SURFACE_FACTOR,
    launch_biodegradation!, launch_photooxidation!, launch_evaporation!,
    N_COMPONENTS
using Zarr

@testset "CICOILPhysics" begin

    @testset "Constants" begin
        @test length(K_BIODEG_ALL)  == 17
        @test length(K_PHOTO_AROM) == 8
        @test all(x -> x > 0, K_BIODEG_ALL)
        @test all(x -> x > 0, K_PHOTO_AROM)
        @test BIODEG_SURFACE_FACTOR ≈ 1.0f0/30.0f0
    end

    @testset "OilProps keyword constructor" begin
        Mw   = fill(200.0f0, 17)
        Tb   = range(350.0f0, 700.0f0, length=17) |> collect .|> Float32
        frac = fill(1.0f0/17.0f0, 17)
        props = OilProps(; Mw, Tb, oil_mass_frac=frac)
        @test props.comp_Mw[1] ≈ 200.0f0
        @test props.Cox_A0[1] ≈ 1.77838f-3 * 350.0f0 + 2.186529f0 atol=1f-4
    end

    @testset "OilParticles allocation" begin
        p = OilParticles(1000)
        @test nparticles(p) == 1000
        @test size(p.mass_components) == (1000, 17)
        @test eltype(p.mass_oil) == Float32
        totals = sum(p.mass_components, dims=2)
        @test all(totals .≈ 185.6f0)
    end

    @testset "Biodegradation — mass conservation" begin
        N = 100
        p = OilParticles(N; initial_mass_kg=1000.0f0)
        fill!(p.z, -5.0f0)
        env = Environment(Dict(
            "sea_water_temperature"   => fill(293.15f0, N),
            "x_sea_water_velocity"    => zeros(Float32, N),
            "y_sea_water_velocity"    => zeros(Float32, N),
            "x_wind"                  => zeros(Float32, N),
            "y_wind"                  => zeros(Float32, N),
            "sea_surface_wave_significant_height" => zeros(Float32, N),
            "sea_surface_wave_period_at_variance_spectral_density_maximum" => zeros(Float32, N),
        ); uv_irradiance=0.0f0)

        initial_mass = sum(p.mass_oil)

        launch_biodegradation!(p, env, 86400.0f0)

        final_oil  = sum(p.mass_oil)
        total_biodeg = sum(p.mass_biodegraded_from_oil)
        @test final_oil + total_biodeg ≈ initial_mass rtol=1f-5
        @test total_biodeg > 0
    end

    @testset "Photooxidation — aromatics only" begin
        N = 50
        p = OilParticles(N; initial_mass_kg=1000.0f0)
        fill!(p.z, 0.0f0)

        env = Environment(Dict(
            "sea_water_temperature"   => fill(302.15f0, N),
            "x_sea_water_velocity"    => zeros(Float32, N),
            "y_sea_water_velocity"    => zeros(Float32, N),
            "x_wind"                  => zeros(Float32, N),
            "y_wind"                  => zeros(Float32, N),
            "sea_surface_wave_significant_height" => zeros(Float32, N),
            "sea_surface_wave_period_at_variance_spectral_density_maximum" => zeros(Float32, N),
        ); uv_irradiance=50.0f0)

        launch_photooxidation!(p, env, 3600.0f0)

        m_phox = sum(p.mass_photooxidized)
        m_op   = sum(p.mass_op_dissolved)
        comps  = p.mass_components

        ali_before = N * (1000.0f0 / 17.0f0) * 8.0f0
        ali_after  = sum(comps[:, 1:8])
        @test ali_after ≈ ali_before rtol=1f-5
        @test m_phox + m_op > 0
        @test m_op ≈ m_phox rtol=0.01f0
    end

    @testset "Evaporation — surface only" begin
        N = 50
        Mw   = fill(200.0f0, 17)
        Tb   = range(350.0f0, 700.0f0, length=17) |> collect .|> Float32
        frac = fill(1.0f0/17.0f0, 17)
        props = OilProps(; Mw, Tb, oil_mass_frac=frac)

        p = OilParticles(N; initial_mass_kg=1000.0f0, oil_props=props)
        fill!(p.z, 0.0f0)
        fill!(p.spillet_thickness, 0.001f0)
        fill!(p.oil_density, 880.0f0)
        fill!(p.oil_molar_mass, 200.0f0)

        env = Environment(Dict(
            "sea_water_temperature"   => fill(303.15f0, N),
            "x_sea_water_velocity"    => zeros(Float32, N),
            "y_sea_water_velocity"    => zeros(Float32, N),
            "x_wind"                  => fill(10.0f0, N),
            "y_wind"                  => zeros(Float32, N),
            "sea_surface_wave_significant_height" => zeros(Float32, N),
            "sea_surface_wave_period_at_variance_spectral_density_maximum" => zeros(Float32, N),
        ); uv_irradiance=0.0f0)

        Kt = compute_transfer_coeff(
            fill(10.0f0, N), zeros(Float32, N),
            fill(200.0f0, N), fill(303.15f0, N),
        )
        initial = sum(p.mass_oil)

        launch_evaporation!(p, env, props, Kt, 3600.0f0)

        evap  = sum(p.mass_evaporated)
        final_m = sum(p.mass_oil)
        @test final_m + evap ≈ initial rtol=1f-5
        @test evap > 0
    end

    @testset "ZarrParticleWriter round-trip" begin
        zarr_path = joinpath(tempdir(), "test_cicoil_particles.zarr")

        N = 20
        n_obs = 3
        p = OilParticles(N; initial_mass_kg=500.0f0)
        fill!(p.z, -5.0f0)

        env = Environment(Dict(
            "sea_water_temperature"   => fill(293.15f0, N),
            "x_sea_water_velocity"    => zeros(Float32, N),
            "y_sea_water_velocity"    => zeros(Float32, N),
            "x_wind"                  => zeros(Float32, N),
            "y_wind"                  => zeros(Float32, N),
            "sea_surface_wave_significant_height" => zeros(Float32, N),
            "sea_surface_wave_period_at_variance_spectral_density_maximum" => zeros(Float32, N),
        ); uv_irradiance=0.0f0)

        writer = ZarrParticleWriter(zarr_path, N, n_obs)

        # Snapshot 1: initial state
        save_snapshot!(writer, p, 0.0, N)

        # Run 1 day of biodegradation
        launch_biodegradation!(p, env, 86400.0f0)

        # Snapshot 2: after biodeg
        save_snapshot!(writer, p, 86400.0, N)

        # Snapshot 3: after another day
        launch_biodegradation!(p, env, 86400.0f0)
        save_snapshot!(writer, p, 172800.0, N)

        # ── Read back and verify ─────────────────────────────────────────────
        g = zopen(Zarr.DirectoryStore(zarr_path))

        # Check shapes
        @test size(g["mass_oil"]) == (N, n_obs)
        @test size(g["mass_components"]) == (N, 17, n_obs)
        @test size(g["obs_time"]) == (n_obs,)

        # Check times
        @test g["obs_time"][1] ≈ 0.0f0
        @test g["obs_time"][2] ≈ 86400.0f0
        @test g["obs_time"][3] ≈ 172800.0f0

        # Check n_active
        @test g["obs_npart"][1] == N
        @test g["obs_npart"][2] == N

        # Mass conservation: initial mass = mass_oil + mass_biodegraded at each snapshot
        for oi in 1:n_obs
            m_oil   = sum(g["mass_oil"][:, oi])
            m_biodeg = sum(g["mass_biodegraded"][:, oi])
            @test m_oil + m_biodeg ≈ N * 500.0f0 rtol=1f-4
        end

        # Biodegradation progressed: snapshot 2 has less oil than snapshot 1
        @test sum(g["mass_oil"][:, 2]) < sum(g["mass_oil"][:, 1])

        # Component field shape
        @test size(g["mass_biodeg_comp"]) == (N, 17, n_obs)

        # Cleanup
        rm(zarr_path, recursive=true)
    end
end
