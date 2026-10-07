# Discrete Prony Series Relaxation Spectrum

using StaticArrays
using LinearAlgebra

"""
    PronySpectrum{N, T<:Real} <: AbstractRelaxationSpectrum

Discrete distribution of relaxation times (Prony series / generalized Maxwell modes):
    Φ(z) = ∑_{i=1}^N g_i * exp( - z / λ_i )
with ∑_{i=1}^N g_i = 1, g_i > 0, λ_i > 0.

# Fields
- `g::SVector{N, T}`: Normalized weighting coefficients for each mode (∑ g_i = 1).
- `lambda::SVector{N, T}`: Dimensionless relaxation time scaling factors for each mode.
"""
struct PronySpectrum{N, T<:Real} <: AbstractRelaxationSpectrum
    g::SVector{N, T}
    lambda::SVector{N, T}

    function PronySpectrum(g::SVector{N, T}, lambda::SVector{N, T}) where {N, T<:Real}
        N > 0 || throw(ArgumentError("Number of modes N must be positive"))
        all(>(zero(T)), lambda) || throw(ArgumentError("All relaxation times λ_i must be strictly positive"))
        all(>=(zero(T)), g) || throw(ArgumentError("All weights g_i must be non-negative"))
        sum_g = sum(g)
        isapprox(sum_g, one(T), atol=1e-5) || throw(ArgumentError("Weights must sum to 1, got sum(g) = $sum_g"))
        # Normalize exactly to 1.0 to eliminate floating point roundoff
        normalized_g = g ./ sum_g
        new{N, T}(normalized_g, lambda)
    end
end

PronySpectrum(g::AbstractVector{<:Real}, lambda::AbstractVector{<:Real}) = begin
    length(g) == length(lambda) || throw(DimensionMismatch("g and lambda must have the same length"))
    N = length(g)
    T = promote_type(eltype(g), eltype(lambda), Float64)
    PronySpectrum(SVector{N, T}(g), SVector{N, T}(lambda))
end

function relaxation_function(s::PronySpectrum{N, T}, z::Real) where {N, T}
    z >= zero(z) || throw(ArgumentError("Material time z must be non-negative, got $z"))
    val = zero(T)
    @inbounds for i in 1:N
        val += s.g[i] * exp(-z / s.lambda[i])
    end
    return val
end
