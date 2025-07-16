module Fminunc

export FminuncStruct, init_fminunc, train

using ..Trainer
using Optim  # Julia's optimization library

"""
    FminuncStruct

Immutable struct storing trainer and optimization parameters.
"""
struct FminuncStruct
    trainer::Trainer.TrainerStruct
    opt_tolerance::Float64
    maxevals::Int
    n_plot::Int
    Xtrain::Matrix{Float64}
    Ytrain::Matrix{Float64}
    opt_options::Optim.Options
end

"""
    init_fminunc(s::Dict{String, Any})

Factory function to initialize FminuncStruct from parameter dictionary.
"""
function init_fminunc(s::Dict{String, Any})
    trainer_obj = Trainer.init_trainer(s)

    opt_tolerance = 1e-10
    maxevals = 5000
    n_plot = 1

    Xtrain = s["Xtrain"]
    Ytrain = s["Ytrain"]

    # Build Optim.Options (approximate translation)
    opt_options = Optim.Options(
        g_tol = opt_tolerance,
        iterations = maxevals * 5,
        f_calls_limit = maxevals,
        show_trace = trainer_obj.is_displayed,
        allow_f_increases = true
    )

    return FminuncStruct(trainer_obj, opt_tolerance, maxevals, n_plot, Xtrain, Ytrain, opt_options)
end

"""
    train(fmin::FminuncStruct)

Runs the optimization using Optim.jl.
Returns updated FminuncStruct with updated trainer containing optimized θ.
"""
function train(fmin::FminuncStruct) # This was not updated in matlab so it probably won't work in the translation
    θ₀ = fmin.trainer.designVariable.thetavec

    # Define objective function compatible with Optim.jl
    function f(θ)
        # Adjust to your framework's cost function call
        return fmin.trainer.objectiveFunction["computeCost"](θ, fmin.Xtrain, fmin.Ytrain) # This should be the problem (inexistant computeCost)
    end

    # Run the optimization (Quasi-Newton by default)
    result = Optim.optimize(f, θ₀, Optim.BFGS(); options = fmin.opt_options)

    θ_opt = Optim.minimizer(result)

    # Update trainer with optimized θ
    updated_trainer = Trainer.update_thetavec(fmin.trainer, θ_opt)

    println("Optimization completed with final cost: $(Optim.minimum(result))")

    # Return updated FminuncStruct with updated trainer
    return FminuncStruct(updated_trainer, fmin.opt_tolerance, fmin.maxevals,
                         fmin.n_plot, fmin.Xtrain, fmin.Ytrain, fmin.opt_options)
end


end