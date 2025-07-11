module Trainer

export AbstractTrainer, TrainerStruct, Create, plotCostRegErr, storeValues!, plotEpsOpt

using Plots
using ..Network.LearnableVariables
using ..CostNN

#abstract type AbstractTrainer end
# Define Trainer structure first
mutable struct TrainerStruct
    objectiveFunction::CostNNStruct
    designVariable::LearnableVars
    xIter::Vector{Vector{Float64}}
    nPlot::Int
    isDisplayed::Bool
    costHist::Matrix{Float64}
    optHist::Matrix{Float64}
    epoch_counter::Int
    #maxEpochs::Int64
end

mutable struct OptInfo
    epsilon::Float64
    gnorm::Float64
end

function TrainerStruct(cParams::Dict{String, Any})
    objFunc = cParams["costFunc"]
    designVar = cParams["designVariable"]
    maxEpochs = cParams["maxEpochs"]

    return TrainerStruct(
        objFunc,
        designVar,
        [zeros(length(designVar.thetavec)) for _ in 1:maxEpochs], # xIter
        1,                              # nPlot
        false,                          # isDisplayed
        zeros(maxEpochs, 3),            # costHist
        zeros(maxEpochs, 2),            # optHist
        1                               # epoch_counter
        #maxEps
    )
end

# Now include trainer submodules
include("SGD/SGD.jl")
include("Fminunc/Fminunc.jl")
include("Nesterov/Nesterov.jl")
include("RMSProp/RMSProp.jl")

# Import their constructors
using .SGD
using .Fminunc
using .Nesterov
using .RMSProp

function plotCostRegErr(t::TrainerStruct, v::Vector{Int})
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

function Create(cParams::Dict{String, Any})
    trainerType = cParams["type"]
    if trainerType == "SGD"
        return SGDStruct(cParams)
    elseif trainerType == "Fminunc"
        return FminuncStruct(cParams)
    elseif trainerType == "Nesterov"
        return NesterovStruct(cParams)
    elseif trainerType == "RMSProp"
        return RMSPropStruct(cParams)
    else
        error("Unknown trainer type: $trainerType")
    end
end

function storeValues!(
    t::TrainerStruct,
    x::Vector{Float64},
    f::Float64,
    opt::OptInfo
)
    # No need for an "if state == init" because t.costHist and t.optHist already initialized to zeros
        epoch = t.epoch_counter

        t.costHist[epoch, :] .= [
            f,
            t.objectiveFunction.regularization,
            t.objectiveFunction.loss
        ]
        t.optHist[epoch, :] .= [opt.gnorm, opt.epsilon]
        t.xIter[epoch] .= x # in-place copy

        t.epoch_counter += 1  # increment counter
end

function plotEpsOpt(t::TrainerStruct, v::Vector{Int})
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