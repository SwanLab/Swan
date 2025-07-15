module RMSProp

export RMSPropStruct, init_RMSProp, step

using ..SGD  # Import SGD module for SGDStruct
using LinearAlgebra

"""
    RMSPropStruct

Immutable struct storing RMSProp optimizer state.
"""
struct RMSPropStruct
    sgd::SGD.SGDStruct
    ρ::Float64
    r::Vector{Float64}
end

"""
    init_RMSProp(s::Dict)

Initializes RMSPropStruct from parameter dictionary.
"""
function init_RMSProp(s::Dict{String, Any})
    sgd_obj = SGD.init_SGD(s)

    ρ = s["rho"]

    θ = s["designVariable"].thetavec
    r = zeros(length(θ))

    return RMSPropStruct(sgd_obj, ρ, r)
end

"""
    step(opt::RMSPropStruct, θ::Vector, ε::Float64, grad::Vector, F, Xb, Yb)

Performs one RMSProp update step. Returns (θ_new, grad_new, updated_optimizer).
"""
function step(opt::RMSPropStruct, θ::Vector{Float64}, ε::Float64, grad::Vector{Float64}, F, Xb, Yb)
    # Update running average of squared gradients
    r_new = opt.ρ .* opt.r .+ (1 .- opt.ρ) .* (grad .^ 2)

    # Parameter update with RMSProp scaling
    θ_new = θ .- ε ./ sqrt.(r_new .+ 1e-6) .* grad

    # Evaluate gradient at new θ if needed (keeping original for standard RMSProp)
    grad_new = grad  # no reevaluation for vanilla RMSProp

    # Return updated optimizer struct
    updated_opt = RMSPropStruct(opt.sgd, opt.ρ, r_new)

    return θ_new, grad_new, updated_opt
end



end