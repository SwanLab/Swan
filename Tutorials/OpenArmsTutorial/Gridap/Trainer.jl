module Trainer

export TrainerStruct, init_Trainer, plot_cost_reg_err, store_values, plot_eps_opt, opt_info

using Plots
using ..Network.LearnableVariables
using ..CostNN

"""
    TrainerStruct

Immutable struct storing training state and history.
"""
struct TrainerStruct
    objectiveFunction::CostNNStruct
    designVariable::LearnableVars
    xIter::Vector{Vector{Float64}}
    nPlot::Int
    isDisplayed::Bool
    costHist::Matrix{Float64}
    optHist::Matrix{Float64}
    epoch_counter::Int
end

"""
    opt_info

Struct storing optimizer info.
"""
struct opt_info
    ε::Float64
    gnorm::Float64
end

"""
    init_Trainer(cParams)

Initializes TrainerStruct from parameter dictionary.
"""
function init_Trainer(cParams::Dict{String, Any})
    objFunc = cParams["costFunc"]
    designVar = cParams["designVariable"]
    maxEpochs = cParams["maxEpochs"]

    # Initialize xIter as empty vectors of correct dimension
    xIter = [zeros(length(designVar.thetavec)) for _ in 1:maxEpochs]

    return TrainerStruct(
        objFunc,
        designVar,
        xIter,
        1,
        false,
        zeros(maxEpochs, 3),
        zeros(maxEpochs, 2),
        1
    )
end

"""
    plot_cost_reg_err(t, v)

Plots cost minimization curve.
"""
function plot_cost_reg_err(t::TrainerStruct, v::Vector{Int})
    idxs = 2:length(v)
    fvals = t.costHist[2:end, 1]

    plot(
        v[idxs], fvals;
        seriestype = :scatter,
        marker = (:diamond, 8, :blue, stroke(0)),
        linestyle = :dash,
        color = :blue,
        xlabel = "Iterations",
        ylabel = "Function Values",
        title = "Cost minimization",
        xlims = (1, maximum(v)),
        legend = false,
        grid = true
    )
end

"""
    store_values(t, thetavec, f, opt)

Returns updated TrainerStruct with stored θ, cost, and optimizer info.
"""
function store_values(t::TrainerStruct, θ::Vector{Float64}, f::Float64, opt::opt_info)
    epoch = t.epoch_counter

    new_costHist = copy(t.costHist)
    new_costHist[epoch, :] .= [f, t.objectiveFunction.regularization, t.objectiveFunction.loss] # Problematic if function is used, no such names

    new_optHist = copy(t.optHist)
    new_optHist[epoch, :] .= [opt.gnorm, opt.ε]

    new_xIter = copy(t.xIter)
    new_xIter[epoch] = copy(θ)

    return TrainerStruct(
        t.objectiveFunction,
        t.designVariable,
        new_xIter,
        t.nPlot,
        t.isDisplayed,
        new_costHist,
        new_optHist,
        epoch + 1
    )
end

"""
    plot_eps_opt(t, v)

Plots step size vs iterations.
"""
function plot_eps_opt(t::TrainerStruct, v::Vector{Int})
    plot(
        v[2:end],
        t.optHist[2:end, 2],
        seriestype = :scatter,
        line = (:dash, :green),
        marker = (:diamond, 6, :green),
        xlabel = "Iterations",
        ylabel = "Learning Rate",
        title = "Step Size vs Iter",
        legend = false,
        xlims = (1, maximum(v))
    )
end


end