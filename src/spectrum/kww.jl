# Kohlrausch-Williams-Watts (KWW) Stretched Exponential Spectrum

"""
    KWWSpectrum{T<:Real} <: AbstractRelaxationSpectrum

Kohlrausch-Williams-Watts (KWW) stretched exponential relaxation response:
    Φ(z) = exp( - z^β )

# Fields
- `beta::T`: Stretching exponent (0 < β ≤ 1). β = 1 corresponds to single exponential Debye relaxation.
"""
struct KWWSpectrum{T<:Real} <: AbstractRelaxationSpectrum
    beta::T

    function KWWSpectrum(beta::Real)
        b = Float64(beta)
        (0.0 < b <= 1.0) || throw(ArgumentError("KWW exponent β must be in (0, 1], got $beta"))
        new{Float64}(b)
    end
end

"""
    relaxation_function(spectrum::AbstractRelaxationSpectrum, z::Real)

Evaluates the normalized relaxation function Φ(z) at dimensionless reduced material time `z ≥ 0`.
"""
function relaxation_function(s::KWWSpectrum, z::Real)
    z >= zero(z) || throw(ArgumentError("Material time z must be non-negative, got $z"))
    iszero(z) && return one(z)
    return exp(-z^(s.beta))
end
