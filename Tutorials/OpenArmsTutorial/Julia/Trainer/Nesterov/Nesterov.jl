module Nesterov

export NesterovStruct, step

using ..SGD
using LinearAlgebra

mutable struct NesterovStruct
    sgd::SGD.SGDStruct
    alpha::Float64
    v::Vector{Float64}
end

#Constructor for NesterovStruct, wrapping an existing SGDStruct and adding Nesterov-specific fields.
function NesterovStruct(s::Dict{String, Any})
    sgd = SGD.SGDStruct(s)
    alpha = s["alpha"] # It's not clear which argument of s it should be, so just rename the right one to "alpha"
    v = zeros(length(sgd.trainer.designVariable["thetavec"]))  # initialize v with correct dimension

    return NesterovStruct(sgd, alpha, v)
end

#Nesterov accelerated gradient descent step.
function step(obj::NesterovStruct, x::Vector{Float64}, e::Float64, grad::Vector{Float64}, F)
    x_hat = x .+ obj.alpha .* obj.v

    # Evaluate gradient at the look-ahead point x_hat
    _, grad_new = F(x_hat)

    # Update velocity
    obj.v = obj.alpha .* obj.v .- e .* grad_new

    # Update parameters
    x = x .+ obj.v

    return x, grad_new
end

end
