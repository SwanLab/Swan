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
        g_tol = optTolerance,           #stop if gradient is small
        iterations = maxevals * 5,      #maximum number of internal iterations
        f_calls_limit = maxevals,       #call limit for the function
        show_trace = trainer_obj.isDisplayed,   #displays progress if enabled
        allow_f_increases = true        #Allows costs to rise at certain times
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

function train(obj::FminuncStruct)
    
    # Retrieve initial weights
    x0 = copy(obj.trainer.designVariable.thetavec)

    # Define objective function AND gradient together
    # Optim.jl's only_fg! interface: F = cost value, G = gradient vector
    function fg!(F, G, theta)
        # Update thetavec with current theta proposed by BFGS
        obj.trainer.designVariable.thetavec .= theta

        # Compute cost and gradient via CostNN
        CostNN.computeFunctionAndGradient(obj.trainer.objectiveFunction, theta)
        
        j  = obj.trainer.objectiveFunction.value
        dj = obj.trainer.objectiveFunction.gradient

        # Fill gradient in-place if requested by Optim.jl
        if G !== nothing
            G .= dj
        end

        # Return cost value if requested by Optim.jl
        if F !== nothing
            return j
        end
    end

    # Run BFGS optimization with correct Optim.jl syntax
    result = Optim.optimize(
        Optim.only_fg!(fg!),
        x0,
        Optim.BFGS(),
        obj.optOptions
    )

    # Check convergence
    if Optim.converged(result)
        println("Optimization converged successfully.")
    else
        println("Warning: optimization did not converge.")
    end

    # Store optimal weights back into the network
    obj.trainer.designVariable.thetavec .= Optim.minimizer(result)

    println("Final cost: $(round(Optim.minimum(result), digits=6))")
    println("Total iterations: $(Optim.iterations(result))")
    println("Total function calls: $(Optim.f_calls(result))")
end

end