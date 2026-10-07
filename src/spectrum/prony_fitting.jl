# Collocation and Fitting of KWW Stretched Exponential to Discrete Prony Series

using LinearAlgebra
using StaticArrays

"""
    kww_to_prony(kww::KWWSpectrum; N::Int=10, z_min::Real=1e-4, z_max::Real=1e4)

Converts a KWW stretched exponential spectrum with exponent `β` into an optimal `N`-mode `PronySpectrum`.
Uses logarithmic collocation across `[z_min, z_max]` solved via Non-Negative Least Squares (NNLS) with unit-sum normalization.
"""
function kww_to_prony(kww::KWWSpectrum; N::Int=10, z_min::Real=1e-4, z_max::Real=1e4)
    beta = kww.beta
    if isapprox(beta, 1.0, atol=1e-6)
        # Single Debye mode
        return PronySpectrum(SVector{1, Float64}(1.0), SVector{1, Float64}(1.0))
    end

    # Logarithmically distributed relaxation times lambda_i
    # When beta is small, the spectrum is wider
    log_min = log10(z_min)
    log_max = log10(z_max)
    log_lambdas = range(log_min, log_max, length=N)
    lambdas = 10.0 .^ log_lambdas

    # Collocation grid (M points, M > N)
    M = max(50, 4 * N)
    log_z = range(log_min - 0.5, log_max + 0.5, length=M)
    z_grid = 10.0 .^ log_z

    # Target KWW values: Φ(z) = exp(-z^β)
    b = exp.(-(z_grid .^ beta))

    # Kernel matrix: A[j, i] = exp(-z_j / λ_i)
    A = zeros(Float64, M, N)
    for j in 1:M
        for i in 1:N
            A[j, i] = exp(-z_grid[j] / lambdas[i])
        end
    end

    # Solve non-negative least squares: min ||A * g - b||^2 s.t. g >= 0
    g = _solve_nnls_projected_gradient(A, b; max_iters=2000, tol=1e-8)
    
    # Normalize weights so sum(g) == 1.0
    sum_g = sum(g)
    if sum_g > 0
        g .= g ./ sum_g
    else
        g .= 1.0 / N
    end

    return PronySpectrum(SVector{N, Float64}(g), SVector{N, Float64}(lambdas))
end

"""
    _solve_nnls_projected_gradient(A, b; max_iters=2000, tol=1e-8)

Solves non-negative least squares `min 0.5 * ||A * x - b||^2 s.t. x >= 0` using accelerated projected gradient descent.
"""
function _solve_nnls_projected_gradient(A::AbstractMatrix{T}, b::AbstractVector{T}; 
                                       max_iters::Int=2000, tol::T=1e-8) where {T<:Real}
    M, N = size(A)
    AtA = A' * A
    Atb = A' * b

    # Estimate Lipschitz constant (largest eigenvalue of AtA) via power iteration
    L = _estimate_max_eigenvalue(AtA)
    step = 1.0 / (L + 1e-12)

    x = ones(T, N) ./ N
    y = copy(x)
    t = one(T)

    for iter in 1:max_iters
        x_old = copy(x)
        # Gradient of 0.5 ||A*y - b||^2 is AtA*y - Atb
        grad = AtA * y - Atb
        # Projected gradient step
        x = max.(zero(T), y - step .* grad)
        
        # Check convergence
        if norm(x - x_old) < tol * (norm(x) + 1e-12)
            break
        end

        # Nesterov momentum update
        t_next = (1 + sqrt(1 + 4 * t^2)) / 2
        gamma = (t - 1) / t_next
        y = x + gamma .* (x - x_old)
        t = t_next
    end

    return x
end

function _estimate_max_eigenvalue(M::AbstractMatrix{T}; iters::Int=30) where {T<:Real}
    n = size(M, 1)
    v = ones(T, n)
    v /= norm(v)
    for _ in 1:iters
        v = M * v
        nrm = norm(v)
        iszero(nrm) && return one(T)
        v /= nrm
    end
    return dot(v, M * v)
end
