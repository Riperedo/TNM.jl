# Adams-Gibbs Configurational Entropy Kinetics

"""
    AdamsGibbsKinetics{T<:Real} <: AbstractRelaxationKinetics

Adams-Gibbs structural relaxation kinetics based on configurational entropy:
    ln(τ(T, Tf)) = ln_A + C / [T * S_c(Tf)]

where S_c(Tf) = S_c0 + ΔCp * ln(Tf / T0).

# Fields
- `ln_A::T`: Pre-exponential logarithmic factor.
- `C::T`: Adams-Gibbs barrier parameter.
- `Sc0::T`: Reference configurational entropy at T0.
- `delta_Cp::T`: Configurational heat capacity increment.
- `T0::T`: Reference temperature.
"""
struct AdamsGibbsKinetics{T<:Real} <: AbstractRelaxationKinetics
    ln_A::T
    C::T
    Sc0::T
    delta_Cp::T
    T0::T

    function AdamsGibbsKinetics(ln_A::Real, C::Real, Sc0::Real, delta_Cp::Real, T0::Real)
        C > 0 || throw(ArgumentError("C must be positive, got $C"))
        Sc0 > 0 || throw(ArgumentError("Sc0 must be positive, got $Sc0"))
        delta_Cp > 0 || throw(ArgumentError("delta_Cp must be positive, got $delta_Cp"))
        T0 > 0 || throw(ArgumentError("T0 must be positive, got $T0"))
        p1, p2, p3, p4, p5 = promote(Float64(ln_A), Float64(C), Float64(Sc0), Float64(delta_Cp), Float64(T0))
        new{typeof(p1)}(p1, p2, p3, p4, p5)
    end
end

AdamsGibbsKinetics(; ln_A::Real, C::Real, Sc0::Real, delta_Cp::Real, T0::Real) = 
    AdamsGibbsKinetics(ln_A, C, Sc0, delta_Cp, T0)

function log_relaxation_time(k::AdamsGibbsKinetics, T::Real, Tf::Real, ::PhysicalUnits)
    T > zero(T) || throw(ArgumentError("Temperature T must be positive"))
    Tf > zero(Tf) || throw(ArgumentError("Fictive temperature Tf must be positive"))
    Sc = k.Sc0 + k.delta_Cp * log(Tf / k.T0)
    Sc > zero(Sc) || throw(ArgumentError("Configurational entropy Sc(Tf) became non-positive ($Sc)"))
    return k.ln_A + k.C / (T * Sc)
end

function log_relaxation_time(k::AdamsGibbsKinetics, theta::Real, theta_f::Real, u::ReducedUnits)
    T = theta * u.T_ref
    Tf = theta_f * u.T_ref
    ln_tau = log_relaxation_time(k, T, Tf, PhysicalUnits())
    return ln_tau - log(u.tau_ref)
end
