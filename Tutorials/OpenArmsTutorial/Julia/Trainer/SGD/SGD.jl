module SGD

export SGDStruct, compute

using ..Trainer  # Use the parent module where TrainerStruct is defined
using Main.CostNN
using Dates  # For timing (like tic/toc)
using Optim
using Plots
using Printf
using LinearAlgebra

mutable struct SGDStruct
    trainer::TrainerStruct

    fvStop::Float64
    lSearchType::String
    MaxEpochs::Int
    maxFunEvals::Int
    optTolerance::Float64
    earlyStop::Int
    timeStop::Float64
    plotter::Any
    svepoch::Int
    fplot::Vector{Float64}
    learningRate::Float64
    elapsedTime::Float64  # total elapsed time in seconds
end

mutable struct KPI
    epoch::Int
    alarm::Float64
    gnorm::Float64
    cost::Float64
end

function SGDStruct(s::Dict{String, Any})
    trainer = Trainer.TrainerStruct(s)

    return SGDStruct(
        trainer,
        1e-4,                             # fvStop
        "static",                         # lSearchType
        s["maxEpochs"],                   # MaxEpochs
        5000,                             # maxFunEvals
        1e-8,                             # optTolerance
        s["maxEpochs"],                   # earlyStop
        Inf,                              # timeStop
        s["plotter"],                     # plotter
        0,                                # svepoch
        fill(0.0, s["maxEpochs"]),        # fplot
        s["learningRate"],                # learningRate
        0.0                               # elapsedTime
    )
end

function compute(obj::SGDStruct)
    start_time = time()

    x0 = obj.trainer.designVariable.thetavec
    optimize(obj, x0)
    
    obj.elapsedTime = time() - start_time
    println("Elapsed time: $(round(obj.elapsedTime, digits=4)) seconds")
end

function plotCostFunc(obj::SGDStruct)
    epoch = 1:length(obj.fplot)
    if !isempty(epoch)
        plt = plot(epoch, obj.fplot;
            xaxis = "Epochs",
            yaxis = "Function Values",
            title = "Cost Function",
            xscale = :log10,
            yscale = :log10,
            linewidth = 1.8,
            legend = false,
            grid = true
        )
        display(plt)
    else
        println("No data to plot in fplot.")
    end
end


function optimize(obj::SGDStruct, th0::Vector{Float64})
    epsilon      = obj.learningRate
    iter         = -1
    funcount     = 0
    alarm        = 0.0
    minTestError = 1.0
    kpi = KPI(1, 0.0, 1.0, 1.0)
    #=
    KPI = Dict(
        :epoch => 1,
        :alarm => 0,
        :gnorm => 1.0,
        :cost  => 1.0
    )
    =#
    theta = th0
    while !isCriteriaMet(obj, kpi)
        state = iter == -1 ? :init : :iter
        newEpoch = true
        moveBatch = true

        while !obj.trainer.objectiveFunction.isBatchDepleted || newEpoch
            f, grad = computeStochasticFunctionAndGradient(obj, theta, moveBatch)
            epsilon, theta, funcount = lineSearch(obj, theta, grad, f, epsilon, funcount)

            funcount += 1
            iter += 1
            newEpoch = false

            kpi.cost  = f
            kpi.gnorm = norm(grad)
            displayIter(obj, iter, funcount, theta, epsilon, state, kpi)
        end

        kpi.epoch += 1
        #kpi_alarm, minTestError = obj.trainer.objectiveFunction["validateES"](alarm, minTestError)
        kpi_alarm, minTestError = CostNN.validateES(obj.trainer.objectiveFunction, alarm, minTestError)
        
        kpi.alarm = kpi_alarm
    end
end

function computeStochasticFunctionAndGradient(obj::SGDStruct, theta::Vector{Float64}, moveBatch::Bool)
    CostNN.setBatchMover!(obj.trainer.objectiveFunction, moveBatch)
    CostNN.computeStochasticFunctionAndGradient!(obj.trainer.objectiveFunction, theta)
    #f    = obj.trainer.objectiveFunction["value"]
    f    = obj.trainer.objectiveFunction.value
    #grad = obj.trainer.objectiveFunction["gradient"]
    grad = obj.trainer.objectiveFunction.gradient
    return f, grad
end

function lineSearch(
    obj::SGDStruct,
    x::Vector{Float64},
    grad::Vector{Float64},
    fOld::Float64,
    e::Float64,
    funcount::Int
)
    moveBatch = false
    F = theta -> computeStochasticFunctionAndGradient(obj, theta, moveBatch)

    type = obj.lSearchType
    xnew = copy(x)  # Default fallback

    if type == "static"
        xnew = step(obj, x, e, grad)

    elseif type == "decay"
        xnew = step(obj, x, e, grad)
        e *= (1 - 1e-3)

    elseif type == "dynamic"
        f = fOld
        while f >= 1.001 * (fOld - e * dot(grad, grad))
            xnew = step(obj, x, e, grad)
            f, _ = F(xnew)
            e /= 2
            funcount += 1
        end
        e *= 5

    elseif type == "fminbnd"
        # Simple scalar minimization over e using Brent's method
        f_e(e1) = F(step(obj, x, e1, grad))[1]
        result = Optim.optimize(f_e, e/10, e*10)
        e = Optim.minimizer(result)
        xnew = step(obj, x, e, grad)

    else
        @warn "Unknown line search type: $type. Defaulting to static step."
        xnew = step(obj, x, e, grad)
    end

    return e, xnew, funcount
end

function updateCriteria(obj::SGDStruct, kpi::KPI)
    return [
        kpi.epoch <= obj.MaxEpochs,
        kpi.alarm < obj.earlyStop,
        kpi.gnorm > obj.optTolerance,
        obj.elapsedTime < obj.timeStop,    # use field now
        kpi.cost > obj.fvStop
    ]
end

function isCriteriaMet(obj::SGDStruct, kpi::KPI)
    criteria = updateCriteria(obj, kpi)
    failedIdx = findfirst(!, criteria)
    itIs = false

    if !isnothing(failedIdx)
        itIs = true
        messages = [
            "Minimization terminated, maximum number of epochs reached $(kpi.epoch)",
            "Minimization terminated, validation set did not decrease in $(obj.earlyStop) epochs",
            "Minimization terminated, reached the optimality tolerance of $(obj.optTolerance)",
            "Minimization terminated, reached the limit time of $(obj.timeStop)",
            "Minimization terminated, reached the target function value of $(obj.fvStop)"
        ]
        if failedIdx ≤ length(messages)
            println(messages[failedIdx])
        end
    end

    return itIs
end

function displayIter(
    obj::SGDStruct,
    iter::Int,
    funcount::Int,
    x::Vector{Float64},
    epsilon::Float64,
    state::Symbol,
    kpi::KPI
)
    opt = Trainer.OptInfo(epsilon * kpi.gnorm, kpi.gnorm)

    printValues(obj, kpi.epoch, funcount, opt, kpi.cost, iter)

    if obj.trainer.isDisplayed && ((kpi.epoch % 25 == 0) || iter == -1)
        obj.trainer.storeValues!(x, kpi.cost, state, opt)
    end
end

function printValues(
    obj::SGDStruct,
    epoch::Int,
    funcount::Int,
    opt::Trainer.OptInfo,
    f::Float64,
    iter::Int
)
    if iter % 20 == 0
        println("                                                        First-order")
        println("Epoch Iteration  Func-count       f(x)        Step-size       optimality")
    end

    @printf("%5d    %5d       %5d    %13.6g  %13.6g   %12.3g\n",
        epoch, iter, funcount, f, opt.epsilon, opt.gnorm)

    if epoch != obj.svepoch
        obj.svepoch = epoch
        #=
        # Ensure fplot has enough capacity
        if length(obj.fplot) < epoch
            resize!(obj.fplot, epoch)
        end
        =#
        obj.fplot[epoch] = f
    end
end



function step(obj::SGDStruct, x::Vector{Float64}, e::Float64, grad::Vector{Float64})
    return x .- e .* grad
end




end