module Sh_Func_L2norm

export ShFuncL2normStruct, init_ShFuncL2norm, compute_stochastic_cost_and_gradient, compute_function_and_gradient

using LinearAlgebra

"""
    ShFuncL2normStruct

Immutable structure storing problem design parameters.
"""
struct ShFuncL2normStruct
    #designVariable::Dict{String, Any}
end

"""
    init_ShFuncL2norm(params)

Initializes a ShFuncL2normStruct from parameter dictionary.
"""
function init_ShFuncL2norm(params::Dict{Any, Any})
    return ShFuncL2normStruct()
end

"""
    compute_stochastic_cost_and_gradient(obj, θ)

Returns (cost, gradient, is_batch_depleted).
For L2 norm, no batching logic is needed so is_batch_depleted is always false.
"""
function compute_stochastic_cost_and_gradient(obj::ShFuncL2normStruct, θ::Vector{Float64}, move_batch=nothing)
    j, dj = compute_function_and_gradient(obj, θ)
    isBD = false
    return j, dj, isBD
end

"""
    compute_function_and_gradient(obj, θ)

Returns (cost, gradient).
"""
function compute_function_and_gradient(obj::ShFuncL2normStruct, θ::Vector{Float64})
    j = compute_cost(θ)
    dj = compute_gradient(θ)
    return j, dj
end

"""
    compute_cost(θ)

Returns 0.5 * ||θ||².
"""
compute_cost(θ::Vector{Float64}) = 0.5 * dot(θ, θ)

"""
    compute_gradient(θ)

Returns θ (gradient of 0.5 * ||θ||²).
"""
compute_gradient(θ::Vector{Float64}) = θ

end
