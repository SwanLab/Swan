module Sh_Func_L2norm

export ShFuncL2norm, computeStochasticCostAndGradient, computeFunctionAndGradient

using LinearAlgebra

mutable struct ShFuncL2norm
    designVariable::Dict{String, Any}
end

function ShFuncL2norm(cParams::Dict{String, Any})
    return ShFuncL2norm(cParams["designVariable"])
end

function computeStochasticCostAndGradient(obj::ShFuncL2norm, x::Vector{Float64}, moveBatch=nothing)
    j, dj = computeFunctionAndGradient(obj, x)
    isBD = false
    return j, dj, isBD
end

function computeFunctionAndGradient(obj::ShFuncL2norm, x::Vector{Float64})
    obj.designVariable["thetavec"] = x
    j = computeCost(obj)
    dj = computeGradient(obj)
    return j, dj
end

function computeCost(obj::ShFuncL2norm)
    theta = obj.designVariable["thetavec"]
    return 0.5 * dot(theta, theta)
    #return 0.5 * theta'*theta   # Equivalent to theta * theta'
end

function computeGradient(obj::ShFuncL2norm)
    return obj.designVariable["thetavec"]
    #return Vector{Float64}(obj.designVariable["thetavec"])
end

end
