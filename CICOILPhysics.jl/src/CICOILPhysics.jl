"""
CICOILPhysics.jl — Weathering kernels for CICOIL v2.

CPU/GPU-agnostic port of the weathering physics from CICOILv2 (OpenDrift
extension, CICESE, K. Kotzakoulakis). When CUDA.jl is loaded, the
CICOILPhysicsCUDAExt extension activates GPU kernels automatically.

The 17-component mass-balance array `mass_components[N, 17]` is the core data
structure: rows are particles, columns are pseudo-components in the order
[ali_G1..ali_G8, arom_G1..arom_G8, residue].
"""
module CICOILPhysics

using Zarr

include("backend.jl")
include("constants.jl")
include("oil_props.jl")
include("particles.jl")
include("kernels/evaporation.jl")
include("kernels/emulsification.jl")
include("kernels/biodegradation.jl")
include("kernels/photooxidation.jl")
include("update.jl")
include("bridge.jl")
include("zarr_output.jl")

export OilProps, OilParticles, Environment, WeatheringConfig
export update_weathering!, compute_transfer_coeff
export nparticles
export ZarrParticleWriter, save_snapshot!
export SimState, py_init_state, py_update_weathering
export AbstractBackend, CPUBackend, backend, set_backend!
export backend_name, has_gpu

end # module
