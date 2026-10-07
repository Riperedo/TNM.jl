# Extended VFTH (Vogel-Fulcher-Tammann-Hesse) Relaxation Kinetics

"""
    VFTHKinetics{T<:Real} <: AbstractRelaxationKinetics

Extended Vogel-Fulcher-Tammann-Hesse (VFTH) relaxation kinetics:
    ln(τ(T, Tf)) = ln_A + B * [x / (T - T_inf) + (1 - x) / (Tf - T_inf)]

# Fields
- `ln_A::T`: Pre-exponential logarithmic factor.
- `B::T`: VFTH activation parameter (in Kelvin).
- `T_inf::T`: Ideal glass transition temperature / Vogel temperature (in Kelvin).
- `x::T`: Non-linearity partitioning parameter (0 < x ≤ 1).
"""
struct VFTHKinetics{T<:Real} <: AbstractRelaxationKinetics
    ln_A::T
    B::T
    T_inf::T
    x::T

    function VFTHKinetics(ln_A::Real, B::Real, T_inf::Real, x::Real)
        B > 0 || throw(ArgumentError("VFTH parameter B must be positive, got $B"))
        (0 < x <= 1) || throw(ArgumentError("Non-linearity parameter x must be in (0, 1], got $x"))
        p1, p2, p3, p4 = promote(Float64(ln_A), Float64(B), Float64(T_inf), Float64(x))
        new{typeof(p1)}(p1, p2, p3, p4)
    end
end

VFTHKinetics(; ln_A::Real, B::Real, T_inf::Real, x::Real) = 
    VFTHKinetics(ln_A, B, T_inf, x)

function log_relaxation_time(k::VFTHKinetics, T::Real, Tf::Real, ::PhysicalUnits)
    (T > k.T_inf) || throw(ArgumentError("Temperature T ($T K) must be strictly above T_inf ($(k.T_inf) K)"))
    (Tf > k.T_inf) || throw(ArgumentError("Fictive temperature Tf ($Tf K) must be strictly above T_inf ($(k.T_inf) K)"))
    return k.ln_A + k.B * (k.x / (T - k.T_inf) + (1 - k.x) / (Tf - k.T_inf))
end

function log_relaxation_time(k::VFTHKinetics, theta::Real, theta_f::Real, u::ReducedUnits)
    theta_inf = k.T_inf / u.T_ref
    (theta > theta_inf) || throw(ArgumentError("Dimensionless T ($theta) must be above θ_inf ($theta_inf)"))
    (theta_f > theta_inf) || throw(ArgumentError("Dimensionless Tf ($theta_f) must be above θ_inf ($theta_inf)"))
    
    B_reduced = k.B / u.T_ref
    ln_tau_tilde = (k.ln_A - log(u.tau_ref)) + B_reduced * (k.x / (theta - theta_inf) + (1 - k.x) / (theta_f - theta_inf))
    return ln_tau_tilde
end
