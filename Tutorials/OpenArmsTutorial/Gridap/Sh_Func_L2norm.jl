module Sh_Func_L2norm

export ShFuncL2norm, 
       computeStochasticCostAndGradient, 
       computeFunctionAndGradient,
       computeCost, 
       computeGradient

using LinearAlgebra

"Stateless structure representing the L2 regularization functional"
struct ShFuncL2norm
    # Placeholder for configuration parameters, if needed in future
end

"Optional constructor accepting a Dict, for compatibility with existing design"
ShFuncL2norm(cParams::Dict{String, Any}) = ShFuncL2norm()

function computeStochasticCostAndGradient(obj::ShFuncL2norm, x::Vector{Float64}, moveBatch=nothing)
    j, dj = computeFunctionAndGradient(obj, x)
    isBD = false
    return j, dj, isBD
end

function computeFunctionAndGradient(obj::ShFuncL2norm, x::Vector{Float64})
    j = computeCost(obj, x)
    dj = computeGradient(obj, x)
    return j, dj
end

function computeCost(::ShFuncL2norm, x::Vector{Float64})
    return 0.5 * dot(x, x)
end

function computeGradient(::ShFuncL2norm, x::Vector{Float64})
    return x
end

end