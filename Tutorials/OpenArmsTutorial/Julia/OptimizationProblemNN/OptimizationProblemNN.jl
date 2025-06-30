module OptimizationProblemNN

export OptimizationProblemNNStruct, solve, getTestData, getNetwork, 
       plotCostFunc, plotBoundary, plotConections, plotConfusionMatrix, plotSurface,
       plotImage, computeError, computeOutputValues, computeGradient

using ..Network
using ..LearnableVariables
using ..LossFunctional
using ..Sh_Func_L2norm
using ..CostNN
using ..Trainer
using ..PlotterNN
using ..Data

mutable struct OptimizationProblemNNStruct
    data::Any
    networkParams::Any
    optimizerParams::Any
    costParams::Any
    network::Any
    costFunc::Any
    designVariable::Any
    optimizer::Any
    plotter::Any
    loss::Any
    regularization::Any
end

function OptimizationProblemNNStruct(cParams::Dict{String,Any})
    obj = new(
        cParams["data"],
        cParams["networkParams"],
        cParams["optimizerParams"],
        cParams["costParams"],
        nothing,    # network placeholder
        nothing,    # costFunc placeholder
        nothing,    # designVariable placeholder
        nothing,    # optimizer placeholder
        nothing,    # plotter placeholder
        nothing,    # loss placeholder
        nothing     # regularization placeholder
    )

    createNetwork(obj)
    createDesignVariable(obj)
    createLossFunctional(obj)
    createRegularizationFunctional(obj)
    createCost(obj)
    createPlotter(obj)
    createOptimizer(obj)

    return obj
end






function createNetwork(obj::OptimizationProblemNNStruct)
    s = deepcopy(obj.networkParams)  # avoid mutating original input
    s["data"] = obj.data
    n = Network.Net(s)
    obj.network = n
end

function createDesignVariable(obj::OptimizationProblemNNStruct)
    dv = Network.getLearnableVariables(obj.network)
    obj.designVariable = dv
end

function createLossFunctional(obj::OptimizationProblemNNStruct)
    s = Dict(
        "network"        => obj.network,
        "designVariable" => obj.designVariable,
        "data"           => obj.data,
        "costType"       => obj.costParams["costType"]
    )
    obj.loss = LossFunctional.LossFunctionalStruct(s)
end



end
