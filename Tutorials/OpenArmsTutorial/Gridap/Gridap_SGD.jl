module SGD

export SGDStruct, init_SGD, compute

using ..Trainer
using ..CostNN
using Dates
using Optim
using Plots
using Printf
using LinearAlgebra

"""
    SGDStruct

Immutable structure storing SGD hyperparameters, trainer, and state history.
"""
struct SGDStruct
    trainer::Trainer.TrainerStruct
    fvStop::Float64
    l_search_type::String
    max_epochs::Int
    max_fun_evals::Int
    opt_tolerance::Float64
    early_stop::Int
    time_stop::Float64
    learning_rate::Float64
    svepoch::Int
    fplot::Vector{Float64}
    elapsed_time::Float64
end

"""
    init_SGD(params)

Initializes SGDStruct from parameter dictionary.
"""
function init_SGD(params::Dict{String, Any})
    trainer = Trainer.init_Trainer(params)

    return SGDStruct(
        trainer,
        1e-4,                       # fvStop
        "static",                   # l_search_type
        params["maxEpochs"],
        5000,                       # max_fun_evals
        1e-8,                       # opt_tolerance
        params["maxEpochs"],
        Inf,                        # time_stop
        params["learningRate"],
        0,                          # svepoch
        fill(0.0, params["maxEpochs"]),  # fplot
        0.0                         # elapsed_time
    )
end

function compute(sgd::SGDStruct, θ::Vector{Float64})
    start_time = time()

    kpi = (epoch=1, alarm=0.0, gnorm=1.0, cost=1.0)

    sgd, _ = _optimize(sgd, θ, kpi, start_time)

    elapsed = time() - start_time
    sgd = update_elapsed_time(sgd, elapsed)

    println("Elapsed time: $(round(elapsed, digits=4)) seconds")
    return sgd
end

function _optimize(sgd::SGDStruct, θ::Vector{Float64}, kpi, start_time::Float64)
    ε = sgd.learning_rate
    min_test_error = 1.0
    iter = -1
    funcount = 0

    while !is_criteria_met(sgd, kpi, time()-start_time)
        new_epoch = true
        move_batch = true
        is_batch_depleted = false
        objFunc = sgd.trainer.objectiveFunction

        while !is_batch_depleted || new_epoch
            # Compute stochastic function and gradient
            f, grad, is_batch_depleted, _= compute_stochastic_function_and_gradient(objFunc, θ, move_batch)

            # Line search
            θ_new, ε, funcount_inc = line_search(sgd, θ, grad, f, ε)
            θ = θ_new

            # Update KPI
            gnorm = norm(grad)
            kpi = (epoch=kpi.epoch, alarm=kpi.alarm, gnorm=gnorm, cost=f)

            # Update fplot and svepoch
            sgd = update_fplot_and_svepoch(sgd, kpi.epoch, f)

            # Display iteration info
            optinfo = Trainer.opt_info(ε * gnorm, gnorm)
            sgd = display_iter(sgd, iter, funcount, θ, kpi, optinfo)

            funcount += 1   # For compute_stochastic_function_and_gradient call
            funcount += funcount_inc   # For line search internal calls
            
            iter += 1
            new_epoch = false
        end

        # Epoch increment
        kpi = (epoch=kpi.epoch+1, alarm=kpi.alarm, gnorm=kpi.gnorm, cost=kpi.cost)

        # Early stopping validation
        alarm, min_test_error = CostNN.validate_ES(objFunc, kpi.alarm, min_test_error, θ)
        kpi = (epoch=kpi.epoch, alarm=alarm, gnorm=kpi.gnorm, cost=kpi.cost)
    end

    return sgd, θ
end

function compute_stochastic_function_and_gradient(objFunc, θ::Vector{Float64}, move_batch::Bool)
    objFunc = CostNN.set_batch_mover(objFunc, move_batch)
    jV, djV, is_batch_depleted, Jc = CostNN.compute_stochastic_function_and_gradient(objFunc, θ)
    return jV, djV, is_batch_depleted, Jc
end

"""
    line_search(sgd, θ, grad, fOld, ε)

Returns updated θ, learning rate, and funcount increment.
"""
function line_search(sgd::SGDStruct, θ::Vector{Float64}, grad::Vector{Float64}, fOld::Float64, ε::Float64)
    type = sgd.l_search_type
    θ_new = copy(θ)
    funcount_inc = 0

    if type == "static"
        θ_new .= step(θ, ε, grad)

    elseif type == "decay"
        θ_new .= step(θ, ε, grad)
        ε *= (1 - 1e-3)

    elseif type == "dynamic"
        f = fOld
        while f >= 1.001 * (fOld - ε * dot(grad, grad))
            θ_new .= step(θ, ε, grad)
            f, _, _, _ = compute_stochastic_function_and_gradient(sgd.trainer.objectiveFunction, θ_new, false)
            ε /= 2
            funcount_inc += 1
        end
        ε *= 5

    elseif type == "fminbnd"
        F = θ -> compute_stochastic_function_and_gradient(sgd.trainer.objectiveFunction, θ, false)[1]
        f_e(e1) = F(step(θ, e1, grad))
        result = Optim.optimize(f_e, ε/10, ε*10)
        ε = Optim.minimizer(result)
        θ_new .= step(θ, ε, grad)

    else
        @warn "Unknown line search type: $type. Defaulting to static step."
        θ_new .= step(θ, ε, grad)
    end

    return θ_new, ε, funcount_inc
end

"""
    is_criteria_met(sgd, kpi, elapsed)

Checks all stopping criteria and prints termination message if met.
"""
function is_criteria_met(sgd::SGDStruct, kpi, elapsed::Float64)
    criteria = [
        kpi.epoch <= sgd.max_epochs,
        kpi.alarm < sgd.early_stop,
        kpi.gnorm > sgd.opt_tolerance,
        elapsed < sgd.time_stop,
        kpi.cost > sgd.fvStop
    ]

    failedIdx = findfirst(!, criteria)
    if isnothing(failedIdx)
        return false
    else
        messages = [
            "Minimization terminated: maximum number of epochs reached ($(kpi.epoch))",
            "Minimization terminated: validation set did not decrease in $(sgd.early_stop) epochs",
            "Minimization terminated: reached the optimality tolerance of $(sgd.opt_tolerance))",
            "Minimization terminated: reached the time limit of $(sgd.time_stop)) seconds",
            "Minimization terminated: reached the target function value of $(sgd.fvStop))"
        ]
        println(messages[failedIdx])
        return true
    end
end

"""
    display_iter(sgd, iter, funcount, θ, kpi, optinfo)

Displays iteration info and stores values if needed.
Returns updated SGDStruct.
"""
function display_iter(sgd::SGDStruct, iter::Int, funcount::Int, θ::Vector{Float64}, kpi, optinfo::Trainer.opt_info)
    if iter % 20 == 0
        println("                                                        First-order")
        println("Epoch Iteration  Func-count       f(x)        Step-size       optimality")
    end

    @printf("%5d    %5d       %5d    %13.6g  %13.6g   %12.3g\n",
        kpi.epoch, iter, funcount, kpi.cost, optinfo.ε, optinfo.gnorm)

    if sgd.trainer.isDisplayed && ((kpi.epoch % 25 == 0) || iter == -1)
        updated_trainer = Trainer.store_values(sgd.trainer, θ, kpi.cost, optinfo)
        sgd = update_trainer(sgd, updated_trainer)
    end

    return sgd
end

"""
    update_trainer(sgd, updated_trainer)

Returns SGDStruct with updated trainer.
"""
update_trainer(sgd::SGDStruct, updated_trainer::Trainer.TrainerStruct) = SGDStruct(updated_trainer, sgd.fvStop, sgd.l_search_type,
    sgd.max_epochs, sgd.max_fun_evals, sgd.opt_tolerance, sgd.early_stop, sgd.time_stop,
    sgd.learning_rate, sgd.svepoch, sgd.fplot, sgd.elapsed_time)

"""
    update_fplot_and_svepoch(sgd, epoch, f)

Returns SGDStruct with updated fplot and svepoch.
"""
function update_fplot_and_svepoch(sgd::SGDStruct, epoch::Int, f::Float64)
    new_fplot = copy(sgd.fplot)
    new_fplot[epoch] = f
    return SGDStruct(sgd.trainer, sgd.fvStop, sgd.l_search_type, sgd.max_epochs,
        sgd.max_fun_evals, sgd.opt_tolerance, sgd.early_stop, sgd.time_stop,
        sgd.learning_rate, epoch, new_fplot, sgd.elapsed_time)
end

"""
    update_elapsed_time(sgd, elapsed)

Returns SGDStruct with updated elapsed_time.
"""
update_elapsed_time(sgd::SGDStruct, elapsed::Float64) = SGDStruct(sgd.trainer, sgd.fvStop, sgd.l_search_type,
    sgd.max_epochs, sgd.max_fun_evals, sgd.opt_tolerance, sgd.early_stop, sgd.time_stop,
    sgd.learning_rate, sgd.svepoch, sgd.fplot, elapsed)

"""
    step(θ, ε, grad)

Returns updated θ after gradient step.
"""
step(θ::Vector{Float64}, ε::Float64, grad::Vector{Float64}) = θ .- ε .* grad

end