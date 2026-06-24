# Backend abstraction for CPU/GPU dispatch.
#
# On CPU: ArrayBackend uses Vector/Matrix and for-loops.
# On GPU: the CICOILPhysicsCUDAExt extension overrides these with CuArray + @cuda.
#
# The extension overloads `make_particles`, `make_env`, and `run_kernel!`.

abstract type AbstractBackend end
struct CPUBackend <: AbstractBackend end

# Global backend — set to CPUBackend by default, overridden by CUDA extension
const ACTIVE_BACKEND = Ref{AbstractBackend}(CPUBackend())

function set_backend!(b::AbstractBackend)
    ACTIVE_BACKEND[] = b
end

backend() = ACTIVE_BACKEND[]

# ── Array allocation dispatch ────────────────────────────────────────────────

alloc_vector(::CPUBackend, ::Type{T}, N::Int) where T = zeros(T, N)
alloc_matrix(::CPUBackend, ::Type{T}, N::Int, M::Int) where T = zeros(T, N, M)
fill_vector(::CPUBackend, val::T, N::Int) where T = fill(val, N)

# Copy from host array to backend array (CPU: no-op copy)
to_backend(::CPUBackend, x::AbstractArray) = Array(x)
to_host(::CPUBackend, x::AbstractArray) = x

# ── Kernel dispatch ──────────────────────────────────────────────────────────
# CPU: simple for-loop over particles
# GPU extension overrides this with @cuda launch

function run_kernel!(::CPUBackend, kernel_body!::Function, N::Int, args...)
    for i in 1:N
        kernel_body!(i, args...)
    end
    return
end

# Synchronize (CPU: no-op; GPU: CUDA.synchronize())
sync_backend!(::CPUBackend) = nothing

# ── Query ────────────────────────────────────────────────────────────────────
has_gpu(::CPUBackend) = false
backend_name(::CPUBackend) = "CPU"
