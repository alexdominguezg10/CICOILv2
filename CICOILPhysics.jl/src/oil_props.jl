# OilProps — per-component physical properties for the 17-component
# pseudo-mixture.  Stored as NTuples (stack-allocated, no heap).
#
# Component order in all arrays: [ali_G1..ali_G8, arom_G1..arom_G8, residue]

struct OilProps
    comp_Mw      ::NTuple{17, Float32}   # g/mol (molecular weight)
    comp_Tb      ::NTuple{17, Float32}   # K     (normal boiling point)
    Cox_A0       ::NTuple{17, Float32}   # Cox vapour-pressure coefficient
    Cox_A1       ::NTuple{17, Float32}
    Cox_A2       ::NTuple{17, Float32}
    oil_mass_frac::NTuple{17, Float32}   # initial mass fraction of each component
end

"""
    OilProps(; Mw, Tb, oil_mass_frac)

Construct OilProps directly from arrays (e.g. when called from the Python bridge
where OpenDrift's FluidProps has already loaded the oil-specific values).
"""
function OilProps(;
    Mw           ::AbstractVector{<:Real},
    Tb           ::AbstractVector{<:Real},
    oil_mass_frac::AbstractVector{<:Real},
)
    length(Mw) == length(Tb) == length(oil_mass_frac) == 17 ||
        throw(ArgumentError("all vectors must have length 17"))

    Tb_f = Float32.(Tb)
    A0 = @. 1.77838f-3  * Tb_f + 2.186529f0
    A1 = @. -1.436577f-7 * Tb_f - 1.908404f-3
    A2 = @. -2.734678f-9 * Tb_f + 2.858042f-6

    OilProps(
        NTuple{17, Float32}(Tuple(Float32.(Mw))),
        NTuple{17, Float32}(Tuple(Tb_f)),
        NTuple{17, Float32}(Tuple(A0)),
        NTuple{17, Float32}(Tuple(A1)),
        NTuple{17, Float32}(Tuple(A2)),
        NTuple{17, Float32}(Tuple(Float32.(oil_mass_frac))),
    )
end
