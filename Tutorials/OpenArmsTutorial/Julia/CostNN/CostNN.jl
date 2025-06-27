module CostNN

export CostNNStruct, computeFunctionAndGradient, computeStochasticFunctionAndGradient,
       obtainNumberFields, getTitleFields, getFields, setBatchMover, validateES

mutable struct CostNNStruct
    value::Float64
    gradient::Vector{Float64}
    isBatchDepleted::Bool

    shapeFunctions::Vector{Any}
    weights::Vector{Float64}
    moveBatch::Bool

    shapeValues::Vector{Any}
end

function CostNNStruct(cParams::Dict{String, Any})
    return CostNNStruct(
        0.0,                             # value
        zeros(0),                        # gradient
        false,                           # isBatchDepleted
        cParams["shapeFunctions"],       # shapeFunctions (vector of modules)
        cParams["weights"],              # weights for each shape function
        true,                            # moveBatch default
        Any[]                            # shapeValues
    )
end

function computeFunctionAndGradient(obj::CostNNStruct, x::Vector{Float64})
    compFunc = (shI, x) -> shI.computeFunctionAndGradient(x)
    jV, djV = computeValueAndGradient(obj, x, compFunc)
    obj.value = jV
    obj.gradient = djV
end




end