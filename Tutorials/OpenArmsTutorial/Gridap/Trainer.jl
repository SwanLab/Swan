module Trainer

export TrainerStruct, init_Trainer, plot_cost_reg_err, store_values, plot_eps_opt, opt_info, update_design_variable

using Plots
using ..Network.LearnableVariables
using ..CostNN

"""
    TrainerStruct

Immutable struct storing training state and history.
"""
struct TrainerStruct
    objective_function::CostNNStruct
    design_variable::LearnableVars
    xIter::Vector{Vector{Float64}}
    nPlot::Int
    is_displayed::Bool
    cost_hist::Matrix{Float64}
    opt_hist::Matrix{Float64}
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
    fvals = t.cost_hist[2:end, 1]

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

    new_cost_hist = copy(t.cost_hist)
    new_cost_hist[epoch, :] .= [f, t.objective_function.regularization, t.objective_function.loss] # Problematic if function is used, no such names

    new_opt_hist = copy(t.opt_hist)
    new_opt_hist[epoch, :] .= [opt.gnorm, opt.ε]

    new_xIter = copy(t.xIter)
    new_xIter[epoch] = copy(θ)

    return TrainerStruct(
        t.objective_function,
        t.design_variable,
        new_xIter,
        t.nPlot,
        t.is_displayed,
        new_cost_hist,
        new_opt_hist,
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
        t.opt_hist[2:end, 2],
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

function update_design_variable(t::TrainerStruct, θ_new::Vector{Float64})
    updated_design_variable = LearnableVariables.update_thetavec(t.design_variable, θ_new)
    
    return TrainerStruct(
        t.objective_function,
        updated_design_variable,
        t.xIter,
        t.nPlot,
        t.is_displayed,
        t.cost_hist,
        t.opt_hist,
        t.epoch_counter
    )
end

end