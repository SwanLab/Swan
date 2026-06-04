module CostNN

export CostNNStruct, init_CostNN, compute_function_and_gradient,
       compute_minibatch_function_and_gradient

using ..LossFunctional
using ..Sh_Func_L2norm
using LinearAlgebra

"""
    CostNNStruct

Struct immutable stockant les fonctions de forme (loss + régularisation) et leurs poids.
shapeFunctions[1] : LossFunctionalStruct  (loss principale)
shapeFunctions[2] : ShFuncL2normStruct    (régularisation L2)
"""
struct CostNNStruct
    shapeFunctions::Tuple{LossFunctional.LossFunctionalStruct, Sh_Func_L2norm.ShFuncL2normStruct}
    weights::Vector{Float64}
end

"""
    init_CostNN(params)

Initialise un CostNNStruct depuis un dictionnaire de paramètres.
"""
function init_CostNN(params::Dict{String, Any})
    sf = params["shapeFunctions"]
    return CostNNStruct(
        (sf[1], sf[2]),
        Float64.(params["weights"])
    )
end

"""
    compute_function_and_gradient(obj, θ)

Coût total et gradient sur le batch complet.
Retourne (J_total, grad_total).
"""
function compute_function_and_gradient(obj::CostNNStruct, θ::Vector{Float64})
    J_loss, g_loss = LossFunctional.compute_function_and_gradient(obj.shapeFunctions[1], θ)
    J_reg,  g_reg  = Sh_Func_L2norm.compute_function_and_gradient(obj.shapeFunctions[2], θ)

    w1, w2 = obj.weights[1], obj.weights[2]
    J_total = w1 * J_loss + w2 * J_reg
    g_total = w1 .* g_loss .+ w2 .* g_reg

    return J_total, g_total
end

"""
    compute_minibatch_function_and_gradient(obj, Xb, Yb, θ)

Coût total et gradient sur un mini-batch (Xb, Yb) déjà extrait.
Retourne (J_total, grad_total).
La responsabilité du découpage en batches appartient à Adam.
"""
# function compute_minibatch_function_and_gradient(obj::CostNNStruct,
#                                                   Xb::Matrix{Float64},
#                                                   Yb::Matrix{Float64},
#                                                   θ::Vector{Float64})
#     J_loss, g_loss = LossFunctional.compute_minibatch_function_and_gradient(
#                          obj.shapeFunctions[1], Xb, Yb, θ)
#     J_reg,  g_reg  = Sh_Func_L2norm.compute_function_and_gradient(obj.shapeFunctions[2], θ)

#     w1, w2 = obj.weights[1], obj.weights[2]
#     J_total = w1 * J_loss + w2 * J_reg
#     g_total = w1 .* g_loss .+ w2 .* g_reg

#     return J_total, g_total
# end

end