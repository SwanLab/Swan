module OptimizationProblemNN

export OptimizationProblemNNStruct, init_OptimizationProblemNN,
       solve, get_test_data, get_network, plot_cost, plot_predictions,
       compute_output_values, compute_gradient

using ..LearnableVariables
using ..Network
using ..LossFunctional
using ..Sh_Func_L2norm
using ..CostNN
using ..Optimizer
using ..PlotterNN
using ..Data

"""
    OptimizationProblemNNStruct

Struct centrale orchestrant toutes les composantes du problème d'optimisation.
"""
struct OptimizationProblemNNStruct
    data::DataStruct
    network::Net
    loss::LossFunctional.LossFunctionalStruct
    reg::Sh_Func_L2norm.ShFuncL2normStruct
    cost::CostNN.CostNNStruct
    optimizer::Union{Optimizer.AdamStruct, Optimizer.SGDStruct}
    plotter::PlotterNN.PlotterNNStruct
end

"""
    init_OptimizationProblemNN(cParams)

Initialise toutes les composantes dans l'ordre de dépendance correct.
"""
function init_OptimizationProblemNN(cParams::Dict{String, Any})
    data             = cParams["data"]
    network_params   = cParams["networkParams"]
    optimizer_params = cParams["optimizerParams"]
    cost_params      = cParams["costParams"]

    # 1. Réseau + variables apprenables
    network        = init_network(merge(network_params, Dict("data" => data)))
    learnable_vars = network.learnable_variables

    # 2. Fonctions de forme
    loss = init_lossfunctional(Dict(
        "network"  => network,
        "data"     => data,
        "costType" => cost_params["costType"]
    ))

    reg = init_ShFuncL2norm(Dict{Any,Any}())

    # 3. Fonction coût totale
    cost = init_CostNN(Dict{String, Any}(
        "shapeFunctions" => [loss, reg],
        "weights"        => [1.0, cost_params["λ"]]
    ))

    # 4. Optimiseur 
    optimizer = _init_optimizer(optimizer_params, cost, learnable_vars)

    # 5. Plotter
    plotter = init_plotter_nn(Dict(
        "data"     => data,
        "network"  => network,
        "costFunc" => cost
    ))

    return OptimizationProblemNNStruct(
        data, network, loss, reg, cost, optimizer, plotter
    )
end

"""
    solve(opt)

Lance l'entraînement. Retourne (optimizer_final, θ_optimal).
"""
function solve(opt::OptimizationProblemNNStruct)
    θ = opt.network.learnable_variables.thetavec
    return Optimizer.compute(opt.optimizer, θ)
end

"""
    get_test_data(opt)

Retourne (Xtest, Ytest).
"""
function get_test_data(opt::OptimizationProblemNNStruct)
    return opt.data.Xtest, opt.data.Ytest
end

"""
    get_network(opt)

Retourne le réseau de neurones.
"""
function get_network(opt::OptimizationProblemNNStruct)
    return opt.network
end

"""
    plot_cost(optimizer)

Affiche la courbe de coût de l'optimiseur Adam.
"""
function plot_cost(optimizer::Union{Optimizer.AdamStruct, Optimizer.SGDStruct})
    Optimizer.plot_cost_func(optimizer)
end

"""
    plot_predictions(opt, θ)

Affiche Ypred vs Ytest sur les données de test.
"""
function plot_predictions(opt::OptimizationProblemNNStruct, θ::Vector{Float64})
    PlotterNN.plot_predictions(opt.plotter, θ)
end

"""
    compute_output_values(opt, X, θ)

Calcule les sorties du réseau pour X avec les paramètres θ.
"""
function compute_output_values(opt::OptimizationProblemNNStruct,
                                X::Matrix{Float64},
                                θ::Vector{Float64})
    return compute_output(opt.network, X, θ)
end

"""
    compute_gradient(opt, X, θ)

Calcule ∂sortie/∂X pour les entrées X.
"""
function compute_gradient(opt::OptimizationProblemNNStruct,
                           X::Matrix{Float64},
                           θ::Vector{Float64})
    return Network.compute_gradient(opt.network, X, θ)
end

function _init_optimizer(optimizer_params::Dict{String, Any},
                          cost::CostNN.CostNNStruct,
                          learnable_vars::LearnableVariables.LearnableVars)
    params = merge(optimizer_params, Dict(
        "costFunc"       => cost,
        "designVariable" => learnable_vars
    ))
    type = optimizer_params["type"]
    if type == "Adam"
        return Optimizer.init_Adam(params)
    elseif type == "SGD"
        return Optimizer.init_SGD(params)
    else
        error("Optimiseur inconnu : $type")
    end
end

end