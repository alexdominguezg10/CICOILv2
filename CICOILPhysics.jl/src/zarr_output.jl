# ZarrParticleWriter — chunked Zarr output for CICOILPhysics particle state.
#
# Pattern from TamocvsOceananigans/plume_eddy_v4.jl: 2D (N, n_obs) arrays
# with Blosc/zstd compression, fill_value=NaN32, per-snapshot callback.

const SCALAR_FIELDS = [
    "z", "diameter", "age_seconds",
    "oil_density", "spillet_thickness",
    "water_fraction", "kinematic_viscosity",
    "mass_oil", "mass_evaporated", "mass_dispersed",
    "mass_biodegraded", "mass_biodegraded_from_oil", "mass_biodegraded_from_water",
    "mass_photooxidized", "mass_op_dissolved", "mass_op_degraded",
    "fraction_evaporated", "fraction_biodegraded_from_oil",
    "fraction_photooxidized", "fraction_op_dissolved", "fraction_op_degraded",
]

const COMPONENT_FIELDS = [
    "mass_components", "mass_evap_comp", "mass_biodeg_comp",
    "mass_photoox_comp", "mass_op_diss_comp", "mass_op_deg_comp",
    "mass_dissolved_comp",
]

const COMPONENT_NAMES = [
    "ali_G1", "ali_G2", "ali_G3", "ali_G4",
    "ali_G5", "ali_G6", "ali_G7", "ali_G8",
    "arom_G1", "arom_G2", "arom_G3", "arom_G4",
    "arom_G5", "arom_G6", "arom_G7", "arom_G8",
    "residue",
]

struct ZarrParticleWriter
    group     ::Zarr.ZGroup{Zarr.DirectoryStore}
    vars      ::Dict{String, Zarr.ZArray}
    comp_vars ::Dict{String, Zarr.ZArray}
    obs_time  ::Zarr.ZArray
    obs_npart ::Zarr.ZArray
    obs_idx   ::Ref{Int}
    n_obs_max ::Int
end

"""
    ZarrParticleWriter(path, N_max, n_obs_max; chunk_particles=min(N_max,1000), chunk_obs=0)

Create a Zarr store for particle weathering output.

- Scalar fields → 2D `(N_max, n_obs_max)` Float32, Blosc/zstd
- Component fields → 3D `(N_max, 17, n_obs_max)` Float32, Blosc/zstd
- `obs_time` / `obs_npart` → 1D `(n_obs_max,)`

`chunk_obs=0` (default) sets chunk_obs = n_obs_max (single chunk in time).
"""
function ZarrParticleWriter(
    path::String,
    N_max::Int,
    n_obs_max::Int;
    chunk_particles::Int = min(N_max, 1000),
    chunk_obs::Int = 0,
)
    chunk_obs = chunk_obs == 0 ? n_obs_max : chunk_obs

    isdir(path) && rm(path, recursive=true)
    pstore = Zarr.DirectoryStore(path)
    pgroup = Zarr.zgroup(pstore; attrs=Dict(
        "scalar_fields"    => SCALAR_FIELDS,
        "component_fields" => COMPONENT_FIELDS,
        "component_names"  => COMPONENT_NAMES,
        "n_components"     => N_COMPONENTS,
        "N_max"            => N_max,
        "n_obs_max"        => n_obs_max,
    ))

    compressor = Zarr.BloscCompressor(cname="zstd", clevel=3)

    vars = Dict{String, Zarr.ZArray}()
    for vname in SCALAR_FIELDS
        vars[vname] = Zarr.zcreate(Float32, pgroup, vname,
            N_max, n_obs_max;
            chunks = (chunk_particles, chunk_obs),
            fill_value = NaN32,
            compressor = compressor)
    end

    comp_vars = Dict{String, Zarr.ZArray}()
    for vname in COMPONENT_FIELDS
        comp_vars[vname] = Zarr.zcreate(Float32, pgroup, vname,
            N_max, N_COMPONENTS, n_obs_max;
            chunks = (chunk_particles, N_COMPONENTS, chunk_obs),
            fill_value = NaN32,
            compressor = compressor)
    end

    obs_time = Zarr.zcreate(Float32, pgroup, "obs_time",
        n_obs_max; fill_value=NaN32)
    obs_npart = Zarr.zcreate(Int32, pgroup, "obs_npart",
        n_obs_max; fill_value=Int32(0))

    ZarrParticleWriter(pgroup, vars, comp_vars, obs_time, obs_npart,
                       Ref(0), n_obs_max)
end

"""
    save_snapshot!(writer, particles, t_seconds, n_active)

Write one observation snapshot of particle state to the Zarr store.
"""
function save_snapshot!(
    w::ZarrParticleWriter,
    p::OilParticles,
    t_seconds::Real,
    n_active::Int,
)
    w.obs_idx[] += 1
    oi = w.obs_idx[]
    oi > w.n_obs_max && return

    w.obs_time[oi]  = Float32(t_seconds)
    w.obs_npart[oi] = Int32(n_active)

    n = min(n_active, nparticles(p))
    n == 0 && return

    # ── Scalar fields ────────────────────────────────────────────────────────
    w.vars["z"][1:n, oi]                          .= @view p.z[1:n]
    w.vars["diameter"][1:n, oi]                    .= @view p.diameter[1:n]
    w.vars["age_seconds"][1:n, oi]                 .= @view p.age_seconds[1:n]
    w.vars["oil_density"][1:n, oi]                 .= @view p.oil_density[1:n]
    w.vars["spillet_thickness"][1:n, oi]            .= @view p.spillet_thickness[1:n]
    w.vars["water_fraction"][1:n, oi]              .= @view p.water_fraction[1:n]
    w.vars["kinematic_viscosity"][1:n, oi]          .= @view p.kinematic_viscosity[1:n]
    w.vars["mass_oil"][1:n, oi]                    .= @view p.mass_oil[1:n]
    w.vars["mass_evaporated"][1:n, oi]             .= @view p.mass_evaporated[1:n]
    w.vars["mass_dispersed"][1:n, oi]              .= @view p.mass_dispersed[1:n]
    w.vars["mass_biodegraded"][1:n, oi]            .= @view p.mass_biodegraded[1:n]
    w.vars["mass_biodegraded_from_oil"][1:n, oi]   .= @view p.mass_biodegraded_from_oil[1:n]
    w.vars["mass_biodegraded_from_water"][1:n, oi] .= @view p.mass_biodegraded_from_water[1:n]
    w.vars["mass_photooxidized"][1:n, oi]          .= @view p.mass_photooxidized[1:n]
    w.vars["mass_op_dissolved"][1:n, oi]           .= @view p.mass_op_dissolved[1:n]
    w.vars["mass_op_degraded"][1:n, oi]            .= @view p.mass_op_degraded[1:n]
    w.vars["fraction_evaporated"][1:n, oi]         .= @view p.fraction_evaporated[1:n]
    w.vars["fraction_biodegraded_from_oil"][1:n, oi].= @view p.fraction_biodegraded_from_oil[1:n]
    w.vars["fraction_photooxidized"][1:n, oi]      .= @view p.fraction_photooxidized[1:n]
    w.vars["fraction_op_dissolved"][1:n, oi]       .= @view p.fraction_op_dissolved[1:n]
    w.vars["fraction_op_degraded"][1:n, oi]        .= @view p.fraction_op_degraded[1:n]

    # ── Component fields (N × 17) ────────────────────────────────────────────
    w.comp_vars["mass_components"][1:n, :, oi]    .= @view p.mass_components[1:n, :]
    w.comp_vars["mass_evap_comp"][1:n, :, oi]     .= @view p.mass_evap_comp[1:n, :]
    w.comp_vars["mass_biodeg_comp"][1:n, :, oi]   .= @view p.mass_biodeg_comp[1:n, :]
    w.comp_vars["mass_photoox_comp"][1:n, :, oi]  .= @view p.mass_photoox_comp[1:n, :]
    w.comp_vars["mass_op_diss_comp"][1:n, :, oi]  .= @view p.mass_op_diss_comp[1:n, :]
    w.comp_vars["mass_op_deg_comp"][1:n, :, oi]   .= @view p.mass_op_deg_comp[1:n, :]
    w.comp_vars["mass_dissolved_comp"][1:n, :, oi] .= @view p.mass_dissolved_comp[1:n, :]

    return
end
