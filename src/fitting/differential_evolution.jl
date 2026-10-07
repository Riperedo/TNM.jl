# Differential Evolution (Genetic Algorithm) Parameter Estimation Engine
# Inspired by Edi-Bot architecture with population tracking and top-K candidate recording

using Random
using LinearAlgebra

"""
    TopKTracker

Tracks and archives the top-K fittest parameter vectors and their loss values during metaheuristic optimization.
"""
mutable struct TopKTracker
    k::Int
    best_params::Vector{Vector{Float64}}
    best_costs::Vector{Float64}

    TopKTracker(k::Int=3) = new(k, Vector{Float64}[], Float64[])
end

function update_tracker!(tracker::TopKTracker, population::Vector{Vector{Float64}}, costs::Vector{Float64})
    # Sort indices by cost in ascending order
    sorted_idx = sortperm(costs)
    top_indices = sorted_idx[1:min(tracker.k, length(sorted_idx))]
    tracker.best_params = [copy(population[i]) for i in top_indices]
    tracker.best_costs = [costs[i] for i in top_indices]
    return tracker
end

"""
    DEResult

Results returned by Differential Evolution parameter optimization.

# Fields
- `best_params::Vector{Float64}`: Globally best parameter vector found.
- `best_cost::Float64`: Cost value of the best parameter vector.
- `top_k_params::Vector{Vector{Float64}}`: Top K parameter vectors.
- `top_k_costs::Vector{Float64}`: Costs of the top K parameter vectors.
- `iterations::Int`: Number of generations executed.
- `nfev::Int`: Total number of objective function evaluations.
- `converged::Bool`: Whether convergence tolerance was met.
"""
struct DEResult
    best_params::Vector{Float64}
    best_cost::Float64
    top_k_params::Vector{Vector{Float64}}
    top_k_costs::Vector{Float64}
    iterations::Int
    nfev::Int
    converged::Bool
end

"""
    optimize_differential_evolution(objective_fn::Function, bounds::Vector{<:Tuple{Real, Real}};
                                    strategy::Symbol=:best1bin,
                                    max_iters::Int=100,
                                    pop_multiplier::Int=15,
                                    F::Union{Real, Tuple{Real, Real}}=(0.5, 1.0),
                                    CR::Real=0.7,
                                    tol::Real=1e-5,
                                    top_k::Int=3,
                                    seed::Union{Nothing, Int}=nothing,
                                    verbose::Bool=false)

Executes Differential Evolution (DE) metaheuristic global optimization to extract TNM model parameters.

# Arguments
- `objective_fn`: Loss function mapping parameter vector `p` -> `loss::Float64`.
- `bounds`: Vector of `(lower_bound, upper_bound)` tuples for each parameter.
- `strategy`: DE mutation strategy (`:best1bin` or `:rand1bin`).
- `max_iters`: Maximum number of generations.
- `pop_multiplier`: Population multiplier (total population size N_p = pop_multiplier * D).
- `F`: Mutation scaling factor (constant Float64 or random interval tuple (F_min, F_max)).
- `CR`: Crossover probability (0 ≤ CR ≤ 1).
- `tol`: Convergence tolerance on population cost spread.
- `top_k`: Number of top candidate parameter sets to track and return.
"""
function optimize_differential_evolution(objective_fn::Function, bounds::Vector{<:Tuple{Real, Real}};
                                        strategy::Symbol=:best1bin,
                                        max_iters::Int=100,
                                        pop_multiplier::Int=15,
                                        F::Union{Real, Tuple{Real, Real}}=(0.5, 1.0),
                                        CR::Real=0.7,
                                        tol::Real=1e-5,
                                        top_k::Int=3,
                                        seed::Union{Nothing, Int}=nothing,
                                        verbose::Bool=false)
    seed !== nothing && Random.seed!(seed)
    
    D = length(bounds)
    N_p = max(4, pop_multiplier * D)
    
    lower_bounds = [Float64(b[1]) for b in bounds]
    upper_bounds = [Float64(b[2]) for b in bounds]
    ranges = upper_bounds .- lower_bounds

    # Initialize population uniformly inside bounds
    population = [lower_bounds .+ rand(D) .* ranges for _ in 1:N_p]
    costs = [objective_fn(population[i]) for i in 1:N_p]
    nfev = N_p

    # Find initial best
    best_idx = argmin(costs)
    gbest_param = copy(population[best_idx])
    gbest_cost = costs[best_idx]

    tracker = TopKTracker(top_k)
    update_tracker!(tracker, population, costs)

    converged = false
    actual_iters = 0

    for iter in 1:max_iters
        actual_iters = iter
        for i in 1:N_p
            # Select mutation factor F
            f_val = if F isa Tuple
                F[1] + rand() * (F[2] - F[1])
            else
                Float64(F)
            end

            # Pick distinct random indices r1, r2, r3 != i
            r1 = rand(1:N_p)
            while r1 == i
                r1 = rand(1:N_p)
            end
            r2 = rand(1:N_p)
            while r2 == i || r2 == r1
                r2 = rand(1:N_p)
            end
            
            # Generate mutant vector v
            v = zeros(Float64, D)
            if strategy === :best1bin
                v = gbest_param .+ f_val .* (population[r1] .- population[r2])
            elseif strategy === :rand1bin
                r3 = rand(1:N_p)
                while r3 == i || r3 == r1 || r3 == r2
                    r3 = rand(1:N_p)
                end
                v = population[r1] .+ f_val .* (population[r2] .- population[r3])
            else
                throw(ArgumentError("Unknown DE strategy: $strategy. Use :best1bin or :rand1bin"))
            end

            # Binomial crossover
            trial = copy(population[i])
            j_rand = rand(1:D) # Ensure at least one parameter is mutated
            for j in 1:D
                if rand() < CR || j == j_rand
                    trial[j] = v[j]
                end
            end

            # Confinement to bounds (clamping)
            trial .= clamp.(trial, lower_bounds, upper_bounds)

            # Evaluate trial vector
            trial_cost = objective_fn(trial)
            nfev += 1

            # Selection
            if trial_cost <= costs[i]
                population[i] = trial
                costs[i] = trial_cost

                if trial_cost < gbest_cost
                    gbest_cost = trial_cost
                    gbest_param = copy(trial)
                end
            end
        end

        update_tracker!(tracker, population, costs)

        if verbose && (iter % 10 == 0 || iter == 1)
            println("DE Generation $iter/$max_iters | Best Cost: $gbest_cost")
        end

        # Check convergence: standard deviation of population costs
        if std(costs) < tol
            converged = true
            break
        end
    end

    return DEResult(gbest_param, gbest_cost, tracker.best_params, tracker.best_costs, actual_iters, nfev, converged)
end
