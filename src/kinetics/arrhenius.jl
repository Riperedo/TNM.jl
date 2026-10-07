# Classical Arrhenius Tool-Narayanaswamy-Moynihan (TNM) Kinetics

"""
    ArrheniusKinetics{T<:Real} <: AbstractRelaxationKinetics

Classical Arrhenius Tool-Narayanaswamy-Moynihan relaxation kinetics:
    ln(τ(T, Tf)) = ln_A + (Δh*/R) * [x/T + (1-x)/Tf]

# Fields
- `ln_A::T`: Natural logarithm of the pre-exponential factor A (where A is in seconds).
- `dh_R::T`: Apparent activation energy divided by universal gas constant Δh*/R (in Kelvin).
- `x::T`: Non-linearity partitioning parameter (0 < x ≤ 1).
"""
struct ArrheniusKinetics{T<:Real} <: AbstractRelaxationKinetics
    ln_A::T
    dh_R::T
    x::T

    function ArrheniusKinetics(ln_A::Real, dh_R::Real, x::Real)
        dh_R > 0 || throw(ArgumentError("dh_R (Δh*/R) must be positive, got $dh_R"))
        (0 < x <= 1) || throw(ArgumentError("Non-linearity parameter x must be in (0, 1], got $x"))
        p1, p2, p3 = promote(Float64(ln_A), Float64(dh_R), Float64(x))
        new{typeof(p1)}(p1, p2, p3)
    end
end

ArrheniusKinetics(; ln_A::Real, dh_R::Real, x::Real) = 
    ArrheniusKinetics(ln_A, dh_R, x)

# Physical units evaluation
function log_relaxation_time(k::ArrheniusKinetics, T::Real, Tf::Real, ::PhysicalUnits)
    T > zero(T) || throw(ArgumentError("Absolute temperature T must be positive, got $T"))
    Tf > zero(Tf) || throw(ArgumentError("Fictive temperature Tf must be positive, got $Tf"))
    return k.ln_A + k.dh_R * (k.x / T + (1 - k.x) / Tf)
end

# Reduced units evaluation
function log_relaxation_time(k::ArrheniusKinetics, theta::Real, theta_f::Real, u::ReducedUnits)
    theta > zero(theta) || throw(ArgumentError("Dimensionless temperature θ must be positive, got $theta"))
    theta_f > zero(theta_f) || throw(ArgumentError("Dimensionless fictive temperature θ_f must be positive, got $theta_f"))
    
    # In reduced units: ε = dh_R / T_ref
    epsilon = k.dh_R / u.T_ref
    ln_tau_tilde = (k.ln_A - log(u.tau_ref)) + epsilon * (k.x / theta + (1 - k.x) / theta_f)
    return ln_tau_tilde
end
