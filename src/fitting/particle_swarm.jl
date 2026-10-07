# Particle Swarm Optimization (PSO) Parameter Estimation Engine
# Inspired by Edi-Bot architecture with swarm bounding and top-K candidate tracking

using Random
using LinearAlgebra

"""
    PSOResult

Results returned by Particle Swarm Optimization.

# Fields
- `best_params::Vector{Float64}`: Globally best parameter vector found.
- `best_cost::Float64`: Cost value of the best parameter vector.
- `top_k_params::Vector{Vector{Float64}}`: Top K parameter sets from personal bests.
- `top_k_costs::Vector{Float64}`: Costs of the top K parameter sets.
- `iterations::Int`: Number of swarm iterations executed.
- `nfev::Int`: Total number of function evaluations.
- `converged::Bool`: Whether convergence tolerance was met.
"""
struct PSOResult
    best_params::Vector{Float64}
    best_cost::Float64
    top_k_params::Vector{Vector{Float64}}
    top_k_costs::Vector{Float64}
    iterations::Int
    nfev::Int
    converged::Bool
end

"""
    optimize_particle_swarm(objective_fn::Function, bounds::Vector{<:Tuple{Real, Real}};
                            max_iters::Int=100,
                            n_particles::Int=40,
                            w::Real=0.5,
                            c1::Real=1.5,
                            c2::Real=1.5,
                            top_k::Int=3,
                            tol::Real=1e-5,
                            seed::Union{Nothing, Int}=nothing,
                            verbose::Bool=false)

Executes Particle Swarm Optimization (PSO) to find optimal TNM model parameters within specified bounds.

# Arguments
- `objective_fn`: Cost evaluation function `p -> Float64`.
- `bounds`: Vector of `(lower, upper)` parameter boundaries.
- `max_iters`: Maximum number of swarm iterations.
- `n_particles`: Number of particles in the swarm.
- `w`: Inertia weight (default: 0.5).
- `c1`: Cognitive acceleration coefficient (default: 1.5).
- `c2`: Social acceleration coefficient (default: 1.5).
- `top_k`: Number of top personal best candidates to archive.
"""
function optimize_particle_swarm(objective_fn::Function, bounds::Vector{<:Tuple{Real, Real}};
                                max_iters::Int=100,
                                n_particles::Int=40,
                                w::Real=0.5,
                                c1::Real=1.5,
                                c2::Real=1.5,
                                top_k::Int=3,
                                tol::Real=1e-5,
                                seed::Union{Nothing, Int}=nothing,
                                verbose::Bool=false)
    seed !== nothing && Random.seed!(seed)

    D = length(bounds)
    lower_bounds = [Float64(b[1]) for b in bounds]
    upper_bounds = [Float64(b[2]) for b in bounds]
    ranges = upper_bounds .- lower_bounds
    v_max = 0.2 .* ranges

    # Initialize positions and velocities
    X = [lower_bounds .+ rand(D) .* ranges for _ in 1:n_particles]
    V = [(-0.1 .* ranges) .+ 0.2 .* rand(D) .* ranges for _ in 1:n_particles]

    pbest = [copy(X[i]) for i in 1:n_particles]
    pbest_costs = [objective_fn(X[i]) for i in 1:n_particles]
    nfev = n_particles

    best_idx = argmin(pbest_costs)
    gbest_pos = copy(pbest[best_idx])
    gbest_cost = pbest_costs[best_idx]

    converged = false
    actual_iters = 0

    for iter in 1:max_iters
        actual_iters = iter
        for i in 1:n_particles
            cost_i = objective_fn(X[i])
            nfev += 1

            if cost_i < pbest_costs[i]
                pbest_costs[i] = cost_i
                pbest[i] = copy(X[i])

                if cost_i < gbest_cost
                    gbest_cost = cost_i
                    gbest_pos = copy(X[i])
                end
            end
        end

        if verbose && (iter % 10 == 0 || iter == 1)
            println("PSO Iteration $iter/$max_iters | Global Best Cost: $gbest_cost")
        end

        # Velocity and position updates
        for i in 1:n_particles
            r1 = rand(D)
            r2 = rand(D)
            
            V[i] .= w .* V[i] .+ c1 .* r1 .* (pbest[i] .- X[i]) .+ c2 .* r2 .* (gbest_pos .- X[i])
            # Velocity clamping
            V[i] .= clamp.(V[i], -v_max, v_max)
            
            X[i] .= X[i] .+ V[i]
            # Position clamping
            X[i] .= clamp.(X[i], lower_bounds, upper_bounds)
        end

        if std(pbest_costs) < tol
            converged = true
            break
        end
    end

    # Extract top K personal bests
    sorted_idx = sortperm(pbest_costs)
    k_sel = min(top_k, length(sorted_idx))
    top_k_params = [copy(pbest[sorted_idx[j]]) for j in 1:k_sel]
    top_k_costs = [pbest_costs[sorted_idx[j]] for j in 1:k_sel]

    return PSOResult(gbest_pos, gbest_cost, top_k_params, top_k_costs, actual_iters, nfev, converged)
end
