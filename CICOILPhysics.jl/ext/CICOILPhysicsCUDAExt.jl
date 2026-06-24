"""
CICOILPhysicsCUDAExt — GPU backend for CICOILPhysics.jl

Loaded automatically when both CICOILPhysics and CUDA are in the environment.
Overrides the backend to use CuArray for particle storage and @cuda for
kernel launches.
"""
module CICOILPhysicsCUDAExt

using CICOILPhysics
using CUDA

import CICOILPhysics: AbstractBackend, alloc_vector, alloc_matrix,
    fill_vector, to_backend, to_host, run_kernel!, sync_backend!,
    has_gpu, backend_name, set_backend!,
    _evaporation_body!, _emulsification_body!,
    _biodegradation_body!, _photooxidation_body!

# ── GPU Backend type ─────────────────────────────────────────────────────────

struct CUDABackend <: AbstractBackend end

alloc_vector(::CUDABackend, ::Type{T}, N::Int) where T = CUDA.zeros(T, N)
alloc_matrix(::CUDABackend, ::Type{T}, N::Int, M::Int) where T = CUDA.zeros(T, N, M)

function fill_vector(::CUDABackend, val::T, N::Int) where T
    v = CUDA.similar(CuArray{T}, N)
    fill!(v, val)
    return v
end

to_backend(::CUDABackend, x::AbstractArray{T}) where T = CuArray{T}(x)
to_host(::CUDABackend, x::CuArray) = Array(x)
to_host(::CUDABackend, x::AbstractArray) = Array(x)

sync_backend!(::CUDABackend) = CUDA.synchronize()
has_gpu(::CUDABackend) = true
backend_name(::CUDABackend) = "CUDA"

# ── CUDA kernel wrappers ────────────────────────────────────────────────────
#
# Each wrapper is a @cuda kernel that calls the _body! function with the
# thread index. The _body! functions use only scalar indexing and basic math,
# so they work on both CPU arrays and CuDeviceArrays.

function _cuda_kernel!(body!, args...)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    body!(i, args...)
    return
end

function run_kernel!(::CUDABackend, kernel_body!::Function, N::Int, args...)
    threads = 256
    blocks = cld(N, threads)
    @cuda blocks=blocks threads=threads _cuda_kernel!(kernel_body!, args...)
    return
end

# ── Auto-activate GPU backend on load ────────────────────────────────────────

function __init__()
    if CUDA.functional()
        set_backend!(CUDABackend())
        @info "CICOILPhysics: CUDA GPU backend activated ($(CUDA.device()))"
    else
        @warn "CICOILPhysics: CUDA loaded but not functional, staying on CPU"
    end
end

end # module
