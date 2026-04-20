module LossFunctional

export LossFunctionalStruct, init_lossfunctional,
       compute_function_and_gradient, compute_minibatch_function_and_gradient

using LinearAlgebra
using ..Network
using ..Data

"""
    LossFunctionalStruct

Conteneur immutable pour les données, le type de coût et le réseau.
"""
struct LossFunctionalStruct
    cost_type::String
    data::DataStruct
    network::Net
end

"""
    init_lossfunctional(params)

Construit un LossFunctionalStruct depuis un dictionnaire de paramètres.
"""
function init_lossfunctional(params::Dict{String, Any})
    return LossFunctionalStruct(
        params["costType"],
        params["data"],
        params["network"]
    )
end

"""
    compute_function_and_gradient(lf, θ)

Coût et gradient sur le batch complet (full-batch).
Retourne (J, grad).
"""
function compute_function_and_gradient(lf::LossFunctionalStruct, θ::Vector{Float64})
    return _forward_and_backward(lf.network, lf.data.Xtrain, lf.data.Ytrain, lf, θ)
end

"""
    compute_minibatch_function_and_gradient(lf, Xb, Yb, θ)

Coût et gradient sur un mini-batch (Xb, Yb) déjà extrait.
Retourne (J, grad).
La responsabilité du découpage en batches appartient à l'optimiseur (Adam).
"""
function compute_minibatch_function_and_gradient(lf::LossFunctionalStruct,
                                                  Xb::Matrix{Float64},
                                                  Yb::Matrix{Float64},
                                                  θ::Vector{Float64})
    return _forward_and_backward(lf.network, Xb, Yb, lf, θ)
end

# ---------------------------------------------------------------------------
# Fonctions internes
# ---------------------------------------------------------------------------

function _forward_and_backward(net::Net, X::Matrix{Float64}, Y::Matrix{Float64},
                                lf::LossFunctionalStruct, θ::Vector{Float64})
    a_vals     = forward_pass(net, X, θ)
    J, dLF     = _compute_cost_and_dLF(lf.cost_type, Y, a_vals[end])
    grad       = backpropagation(net, Y, dLF, a_vals, θ)
    return J, grad
end

"""
    _compute_cost_and_dLF(cost_type, Y, Ŷ)

Retourne (J, dJ/dŶ) selon le type de coût.
"""
function _compute_cost_and_dLF(cost_type::String, Y::Matrix{Float64}, Ŷ::Matrix{Float64})
    if cost_type == "L2"
        diff = Ŷ .- Y
        J    = sqrt(sum(diff .^ 2))
        dLF  = diff
        return J, dLF

    elseif cost_type == "-loglikelihood"
        yp  = Ŷ .- 1e-11   # stabilité numérique uniquement ici
        c   = (1.0 .- Y) .* (-log.(1.0 .- yp)) .+ Y .* (-log.(yp))
        J   = mean(sum(c, dims=2))
        dLF = (yp .- Y) ./ (yp .* (1.0 .- yp))
        return J, dLF

    else
        error("Type de coût invalide : $cost_type")
    end
end

end