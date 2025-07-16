module OptimizationProblemNN

export OptimizationProblemNNStruct, solve, getTestData, getNetwork, 
       plotCostFunc, plotBoundary, plotConections, plotConfusionMatrix, plotSurface,
       plotImage, computeError, computeOutputValues, computeGradient

using ..Network
using ..LossFunctional
using ..Sh_Func_L2norm
using ..CostNN
using ..Trainer
using ..SGD
using ..Nesterov
using ..RMSProp
using ..Fminunc
using ..PlotterNN
using ..Data
# Performance optimization
using Profile 
using ProfileView
using InteractiveUtils

mutable struct OptimizationProblemNNStruct
    data::Data.DataStruct
    networkParams::Dict{String, Any}
    optimizerParams::Dict{String, Any}
    costParams::Dict{String, Any}
    network::Union{Network.Net, Nothing}
    costFunc::Union{CostNN.CostNNStruct, Nothing}
    designVariable::Union{Network.LearnableVariables.LearnableVars, Nothing}
    #optimizer::Union{Trainer.SGD.SGDStruct, Nothing}
    optimizer::Any
    plotter::Union{PlotterNN.PlotterNNStruct, Nothing}
    loss::Union{LossFunctional.LossFunctionalStruct, Nothing}
    regularization::Union{Sh_Func_L2norm.ShFuncL2norm, Nothing}
end

function OptimizationProblemNNStruct(cParams::Dict{String,Any})
    obj = OptimizationProblemNNStruct(
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

# Public methods

function solve(obj::OptimizationProblemNNStruct)
    #@code_warntype 
    Trainer.SGD.compute(obj.optimizer)
end

function getTestData(obj::OptimizationProblemNNStruct)
    #return obj.data["Xtest"], obj.data["Ytest"]
    return obj.data.Xtest, obj.data.Ytest
end

function getNetwork(obj::OptimizationProblemNNStruct)
    return obj.network
end

function plotCostFnc(obj::OptimizationProblemNNStruct)
    Trainer.SGD.plotCostFunc(obj.optimizer)
end

function plotBoundary(obj::OptimizationProblemNNStruct, type::String="filledS") # !!! Cannot work without a method Data.buildModel()
    PlotterNN.plotBoundary(obj.plotter, type)
end

function plotConections(obj::OptimizationProblemNNStruct)
    PlotterNN.plotNetworkStatus(obj.plotter)
end

function plotConfusionMatrix(obj::OptimizationProblemNNStruct)
    PlotterNN.drawConfusionMat(obj.plotter)
end

function plotSurface(obj::OptimizationProblemNNStruct)
    PlotterNN.drawSurfaceResults(obj.plotter)
end

function plotImage(obj::OptimizationProblemNNStruct, row::Int)
    PlotterNN.image(obj.plotter, row)
end

function computeError(obj::OptimizationProblemNNStruct)
    return LossFunctional.getTestError(obj.loss) # Different to the Matlab version since forwardprop() does not exist anymore
end

function computeOutputValues(obj::OptimizationProblemNNStruct, X::Matrix{Float64})
    return Network.computeYOut(obj.network, X)
end

function computeGradient(obj::OptimizationProblemNNStruct, X::Matrix{Float64})
    return Network.networkGradient(obj.network, X)
end

# Private methods

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

function createRegularizationFunctional(obj::OptimizationProblemNNStruct)
    s = Dict("designVariable" => obj.designVariable)
    obj.regularization = Sh_Func_L2norm.ShFuncL2norm(s)
end

function createCost(obj::OptimizationProblemNNStruct)
    s = Dict{String, Any}(
        "shapeFunctions" => [obj.loss, obj.regularization],
        "weights" => [1.0, obj.costParams["lambda"]]
    )
    obj.costFunc = CostNN.CostNNStruct(s)
end

function createOptimizer(obj::OptimizationProblemNNStruct)
    s = Dict{String, Any}(copy(obj.optimizerParams))

    s["costFunc"] = obj.costFunc
    s["designVariable"] = Network.getLearnableVariables(obj.network)
    s["type"] = "SGD" # Specify your optimizer here
    s["data"] = obj.data
    s["Xtrain"] = obj.data.Xtrain
    s["Ytrain"] = obj.data.Ytrain
    s["Xtest"] = obj.data.Xtest
    s["Ytest"] = obj.data.Ytest
    s["plotter"] = obj.plotter
    obj.optimizer = Trainer.Create(s)
end

function createPlotter(obj::OptimizationProblemNNStruct)
    s = Dict(
        "network" => obj.network,
        "data" => obj.data,
        "costFunc" => obj.costFunc
    )
    obj.plotter = PlotterNN.PlotterNNStruct(s)
    println(obj.network.neuronsPerLayer)
end

end
