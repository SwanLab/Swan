module Nesterov

export NesterovStruct, init_nesterov, step

using ..SGD
using LinearAlgebra

"""
    NesterovStruct

Immutable struct storing SGD parameters plus Nesterov-specific parameters.
- sgd: SGDStruct instance
- alpha: momentum coefficient
- v: velocity vector
"""
struct NesterovStruct
    sgd::SGD.SGDStruct
    alpha::Float64
    v::Vector{Float64}
end

"""
    init_nesterov(s::Dict{String,Any})

Factory function to initialize a NesterovStruct from a parameter dictionary.
"""
function init_nesterov(s::Dict{String, Any})
    sgd = SGD.init_sgd(s)
    alpha = s["alpha"]
    v = zeros(length(sgd.trainer.design_variable.thetavec))
    return NesterovStruct(sgd, alpha, v)
end

"""
    step(obj::NesterovStruct, x::Vector{Float64}, e::Float64, grad::Vector{Float64}, F)

Performs one Nesterov accelerated gradient descent step.
Returns updated x, grad_new, and updated NesterovStruct with new velocity.
"""
function step(obj::NesterovStruct, θ::Vector{Float64}, ε::Float64, grad::Vector{Float64}, F)
    θ_hat = θ .+ obj.alpha .* obj.v

    # Evaluate gradient at the look-ahead point x_hat
    _, grad_new = F(θ_hat)

    # Update velocity (no mutation)
    v_new = obj.alpha .* obj.v .- ε .* grad_new

    # Update parameters
    θ_new = θ .+ v_new

    # Return updated x, grad, and updated struct
    return θ_new, grad_new, NesterovStruct(obj.sgd, obj.alpha, v_new)
end

end