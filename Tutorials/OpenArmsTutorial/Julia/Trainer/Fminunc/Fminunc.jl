module Fminunc

export FminuncStruct, train

using ..Trainer
using Optim  # Julia's optimization library

mutable struct FminuncStruct
    trainer::Trainer.TrainerStruct
    optTolerance::Float64
    maxevals::Int
    nPlot::Int
    Xtrain::Matrix{Float64}
    Ytrain::Matrix{Float64}
    optOptions::Optim.Options
end

#Constructor replicating MATLAB's Fminunc initialization.
function FminuncStruct(s::Dict{String, Any})
    trainer_obj = Trainer.TrainerStruct(s)

    optTolerance = 1e-10
    maxevals = 5000
    nPlot = 1

    Xtrain = s["Xtrain"]
    Ytrain = s["Ytrain"]

    # Build Optim.Options (approximate translation)
    optOptions = Optim.Options(
        g_tol = optTolerance,
        iterations = maxevals * 5,
        f_calls_limit = maxevals,
        show_trace = trainer_obj.isDisplayed,
        allow_f_increases = true
    )

    return FminuncStruct(trainer_obj, optTolerance, maxevals, nPlot, Xtrain, Ytrain, optOptions)
end

#Runs the optimization using Optim.jl.
function train(obj::FminuncStruct)
    x0 = obj.trainer.designVariable.thetavec

    # Define objective function compatible with Optim.jl
    function f(theta)
        # Assuming computeCost is accessible; adjust to your framework call
        return obj.trainer.objectiveFunction["computeCost"](theta, obj.Xtrain, obj.Ytrain)
    end

    # Run the optimization (Quasi-Newton by default in Optim.jl)
    result = Optim.optimize(f, x0, Optim.BFGS(); options = obj.optOptions)

    # Store the result if needed
    obj.trainer.designVariable.thetavec = Optim.minimizer(result)

    println("Optimization completed with final cost: $(Optim.minimum(result))")
end

end