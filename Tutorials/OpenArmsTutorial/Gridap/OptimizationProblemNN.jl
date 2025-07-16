module OptimizationProblemNN

export OptimizationProblemNNStruct, init_OptimizationProblemNN,
       solve, get_test_data, get_network, plot_cost,
       plot_boundary, plot_connections, plot_confusion_matrix,
       plot_surface, plot_image, compute_error,
       compute_output_values, compute_gradient

using ..LearnableVariables
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

struct OptimizationProblemNNStruct
    data::DataStruct
    network_params::Dict{String, Any}
    optimizer_params::Dict{String, Any}
    cost_params::Dict{String, Any}
    network::Net
    learnable_vars::LearnableVariables.LearnableVars
    loss::LossFunctional.LossFunctionalStruct
    reg::Sh_Func_L2norm.ShFuncL2normStruct
    cost::CostNN.CostNNStruct
    trainer::Trainer.TrainerStruct
    optimizer::Any
    plotter::PlotterNN.PlotterNNStruct
end

function init_OptimizationProblemNN(cParams::Dict{String, Any})
    data = cParams["data"]
    network_params = cParams["networkParams"]
    optimizer_params = cParams["optimizerParams"]
    cost_params = cParams["costParams"]

    network = init_network(merge(deepcopy(network_params), Dict("data" => data)))

    learnable_vars = get_learnable_variables(network)

    loss = init_lossfunctional(Dict(
        "network" => network,
        "designVariable" => learnable_vars,
        "data" => data,
        "costType" => cost_params["costType"]
    ))

    reg = init_ShFuncL2norm(Dict())

    cost = init_CostNN(Dict{String, Any}(
        "shapeFunctions" => [loss, reg],
        "weights" => [1.0, cost_params["λ"]]
    ))

    trainer = init_Trainer(Dict(
        "costFunc" => cost,
        "designVariable" => learnable_vars,
        "maxEpochs" => optimizer_params["maxEpochs"]
    ))

    optimizer = _init_optimizer(optimizer_params, cost, learnable_vars)

    plotter = init_plotter_nn(Dict(
        "data" => data,
        "network" => network,
        "costFunc" => cost
    ))

    return OptimizationProblemNNStruct(
        data,
        network_params,
        optimizer_params,
        cost_params,
        network,
        learnable_vars,
        loss,
        reg,
        cost,
        trainer,
        optimizer,
        plotter
    )
end
#=
function solve(opt::OptimizationProblemNNStruct)
    compute(opt.optimizer, opt.learnable_vars)
end
=#

function solve(opt::OptimizationProblemNNStruct)
    θvec = opt.learnable_vars.thetavec
    optimizer, θ = compute(opt.optimizer, θvec)
    return optimizer, θ
end

function get_test_data(opt::OptimizationProblemNNStruct)
    return opt.data.Xtest, opt.data.Ytest
end

function get_network(opt::OptimizationProblemNNStruct)
    return opt.network
end

function plot_cost(optimizer::Any)
    plot_cost_func(optimizer)
end

function plot_boundary(opt::OptimizationProblemNNStruct, type::String="filledS")
    PlotterNN.plot_boundary(opt.plotter, type)
end

function plot_connections(opt::OptimizationProblemNNStruct)
    PlotterNN.plot_network_status(opt.plotter)
end

function plot_confusion_matrix(opt::OptimizationProblemNNStruct)
    PlotterNN.draw_confusion_mat(opt.plotter)
end

function plot_surface(opt::OptimizationProblemNNStruct)
    PlotterNN.draw_surface_results(opt.plotter)
end

function plot_image(opt::OptimizationProblemNNStruct, row::Int)
    PlotterNN.image(opt.plotter, row)
end

function compute_error(opt::OptimizationProblemNNStruct, θ::Vector{Float64})
    return LossFunctional.get_test_error(opt.loss, θ)
end

function compute_output_values(opt::OptimizationProblemNNStruct, X::Matrix{Float64}, θ::Vector{Float64})
    return compute_output(opt.network, X, θ)
end

function compute_gradient(opt::OptimizationProblemNNStruct, X::Matrix{Float64}, θ::Vector{Float64})
    return network_gradient(opt.network, X, θ)
end

function _init_optimizer(optimizer_params::Dict{String, Any}, cost::CostNNStruct, learnable_vars::LearnableVars)
    optimizer_type = optimizer_params["type"]
    params = merge(deepcopy(optimizer_params), Dict(
        "costFunc" => cost,
        "designVariable" => learnable_vars
    ))

    if optimizer_type == "SGD"
        return init_SGD(params)
    elseif optimizer_type == "Nesterov"
        return init_Nesterov(params)
    elseif optimizer_type == "RMSProp"
        return init_RMSProp(params)
    elseif optimizer_type == "Fminunc"
        return init_Fminunc(params)
    else
        error("Unknown optimizer type: $optimizer_type")
    end
end

end