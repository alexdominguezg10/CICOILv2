# setup_envs.jl — Initialize CPU and GPU Julia environments
#
# Run once after cloning:
#   julia envs/setup_envs.jl
#
# This creates separate Manifest.toml files for each environment
# so precompile caches never conflict between CPU and GPU nodes.

using Pkg

root = dirname(dirname(@__FILE__))
pkg_path = root  # CICOILPhysics.jl source

for env in ["cpu", "gpu"]
    env_path = joinpath(root, "envs", env)
    println("\n══════════════════════════════════════")
    println("Setting up: envs/$env")
    println("══════════════════════════════════════")

    Pkg.activate(env_path)
    Pkg.develop(path=pkg_path)
    Pkg.instantiate()
    Pkg.precompile()

    println("✓ envs/$env ready")
end

println("\n══════════════════════════════════════")
println("Both environments ready.")
println("")
println("Usage:")
println("  GPU node:  julia --project=envs/gpu")
println("  CPU node:  julia --project=envs/cpu")
println("══════════════════════════════════════")
