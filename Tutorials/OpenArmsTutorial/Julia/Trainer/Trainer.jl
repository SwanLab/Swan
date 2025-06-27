module Trainer

export AbstractTrainer, TrainerStruct, Create, plotCostRegErr, storeValues!, plotEpsOpt

using Plots

# Include trainer submodules
include("SGD/SGD.jl")
include("Fminunc/Fminunc.jl")
include("Nesterov/Nesterov.jl")
include("RMSProp/RMSProp.jl")

# Import their constructors
using .SGD
using .Fminunc
using .Nesterov
using .RMSProp

#abstract type AbstractTrainer end

mutable struct TrainerStruct
    objectiveFunction::Any
    designVariable::Any
    xIter::Vector{Vector{Float64}}
    nPlot::Int
    isDisplayed::Bool
    costHist::Matrix{Float64}
    optHist::Matrix{Float64}
end

function TrainerStruct(cParams::Dict{String, Any})
    objFunc = cParams["costFunc"]
    designVar = cParams["designVariable"]

    return TrainerStruct(
        objFunc,
        designVar,
        Vector{Vector{Float64}}(),  # xIter
        1,                          # nPlot
        false,                      # isDisplayed
        zeros(0, 3),                # costHist
        zeros(0, 2)                 # optHist
    )
end

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
    state::Symbol,
    opt::Dict{Symbol, Float64}
)
    if state == :init
        t.costHist = [zeros(3)]      # equivalent to [0,0,0]
        t.optHist  = [zeros(2)]      # equivalent to [0,0]
        # Julia doesn't need figure handles for plotting
    elseif state == :iter
        cV = [
            f,
            t.objectiveFunction.regularization,
            t.objectiveFunction.loss
        ]
        t.costHist = vcat(t.costHist, cV')
        t.optHist  = vcat(t.optHist, oV')
        oV = [opt[:gnorm], opt[:epsilon]]
        push!(t.optHist, oV)
        

        push!(t.xIter, copy(x))
    else
        @warn "Unknown state passed to storeValues: $state"
    end
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