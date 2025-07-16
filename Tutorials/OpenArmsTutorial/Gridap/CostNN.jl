module CostNN

export CostNNStruct, init_CostNN, compute_function_and_gradient, compute_stochastic_function_and_gradient,
       obtain_number_fields, get_title_fields, get_fields, set_batch_mover, validate_ES

using ..LossFunctional
using ..Sh_Func_L2norm
using ..Network
using LinearAlgebra
using Random

"""
    CostNNStruct

Immutable struct storing shape functions and weights.
"""
struct CostNNStruct
    shapeFunctions::Vector{Any}
    weights::Vector{Float64}
    moveBatch::Bool
end

"""
    init_CostNN(params)

Initializes CostNNStruct from parameter dictionary.
"""
function init_CostNN(params::Dict{String, Any})
    shapeFuncs = params["shapeFunctions"]
    weights = params["weights"]
    return CostNNStruct(shapeFuncs, weights, true)
end

function compute_function_and_gradient(obj::CostNNStruct, θ::Vector{Float64})
    nF = length(obj.shapeFunctions)
    Jc = Vector{Any}(undef, nF)
    dJc = Vector{Any}(undef, nF)

    for iF in 1:nF
        shI = obj.shapeFunctions[iF]
        if shI isa LossFunctional.LossFunctionalStruct
            Jc[iF], dJc[iF] = LossFunctional.compute_function_and_gradient(shI, θ)
        elseif shI isa Sh_Func_L2norm.ShFuncL2normStruct
            Jc[iF], dJc[iF] = Sh_Func_L2norm.compute_function_and_gradient(shI, θ)
        else
            error("Unsupported shape function type: $(typeof(shI))")
        end
    end

    jV = sum(obj.weights .* map(x -> x, Jc))
    djV = sum([obj.weights[i] * dJc[i] for i in 1:nF])

    return jV, djV, Jc
end

function compute_stochastic_function_and_gradient(obj::CostNNStruct, θ::Vector{Float64})
    nF = length(obj.shapeFunctions)
    Jc = Vector{Any}(undef, nF)
    dJc = Vector{Any}(undef, nF)
    bDa = falses(nF)
    for iF in 1:nF
        shI = obj.shapeFunctions[iF]
        if shI isa LossFunctional.LossFunctionalStruct
            Jc[iF], dJc[iF], bDa[iF] = LossFunctional.compute_stochastic_cost_and_gradient(shI, θ, 1, randperm(size(shI.data.Xtrain, 1)), obj.moveBatch)
        elseif shI isa Sh_Func_L2norm.ShFuncL2normStruct
            Jc[iF], dJc[iF], bDa[iF] = Sh_Func_L2norm.compute_stochastic_cost_and_gradient(shI, θ, obj.moveBatch)
        else
            error("Unsupported shape function type: $(typeof(shI))")
        end
    end
    
    jV = sum(obj.weights .* map(x -> x, Jc))
    djV = sum([obj.weights[i] * dJc[i] for i in 1:nF])

    isBatchDepleted = any(bDa)

    return jV, djV, isBatchDepleted, Jc
end

"""
    obtain_number_fields(obj)

Returns number of shape functions.
"""
obtain_number_fields(obj::CostNNStruct) = length(obj.shapeFunctions)

"""
    get_title_fields(obj)

Returns a vector of titles for each shape function with weights.
"""
function get_title_fields(obj::CostNNStruct)
    nF = length(obj.shapeFunctions)
    titles = Vector{String}(undef, nF)

    for iF in 1:nF
        wI = obj.weights[iF]
        titles[iF] = "ShapeFunction $iF (w=$wI)"
    end
    return titles
end

"""
    get_fields(Jc, i)

Returns the i-th shape value from precomputed Jc.
"""
get_fields(Jc::Vector, i::Int) = Jc[i]

"""
    set_batch_mover(obj, moveBatch)

Returns a new CostNNStruct with updated moveBatch flag.
"""
function set_batch_mover(obj::CostNNStruct, moveBatch::Bool)
    return CostNNStruct(obj.shapeFunctions, obj.weights, moveBatch)
end

"""
    validate_ES(obj, alarm, minTestError, θ)

Checks test error and updates alarm/minTestError accordingly.
"""
function validate_ES(obj::CostNNStruct, alarm::Float64, minTestError::Float64, θ::Vector{Float64})
    shI = obj.shapeFunctions[1]

    testError = shI isa LossFunctional.LossFunctionalStruct ? LossFunctional.get_test_error(shI, θ) :
                shI isa Sh_Func_L2norm.ShFuncL2normStruct ? 0.0 : # No test error for L2norm
                error("Unsupported shape function type: $(typeof(shI))")

    if testError < minTestError
        minTestError = testError
        alarm = 0.0
    elseif testError == minTestError
        alarm += 0.5
    else
        alarm += 1.0
    end
    return alarm, minTestError
end


end