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

function computeStochasticFunctionAndGradient!(obj::CostNNStruct, x::Vector{Float64})
    compFunc = (shI, x, moveBatch) -> shI.computeStochasticCostAndGradient(shI, x, moveBatch)
    jV, djV = computeValueAndGradient(obj, x, compFunc)
    obj.value = jV
    obj.gradient = djV
end


function obtainNumberFields(obj::CostNNStruct)
    return length(obj.shapeFunctions)
end

function getTitleFields(obj::CostNNStruct)
    nF = length(obj.shapeFunctions)
    titles = Vector{String}(undef, nF)
    for iF in 1:nF
        wI = obj.weights[iF]
        titleF = obj.shapeFunctions[iF].getTitleToPlot()
        titles[iF] = "$titleF (w=$wI)"
    end
    return titles
end

function getFields(obj::CostNNStruct, i::Int)
    return obj.shapeValues[i]
end

function setBatchMover!(obj::CostNNStruct, moveBatch::Bool)
    obj.moveBatch = moveBatch
end

function validateES(obj::CostNNStruct, alarm::Float64, minTestError::Float64)
    testError = obj.shapeFunctions[1].getTestError(obj.shapeFunctions[1])
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



function computeValueAndGradient(obj::CostNNStruct, x::Vector{Float64}, compFunc)
    nF = length(obj.shapeFunctions)
    bDa = falses(nF)
    Jc = Vector{Any}(undef, nF)
    dJc = Vector{Any}(undef, nF)

    for iF in 1:nF
        shI = obj.shapeFunctions[iF]

        # Try calling with 3 args (stochastic), fallback to 2 args (full batch)
        try
            j, dJ, bD = compFunc(shI, x, obj.moveBatch)
            bDa[iF] = bD
        catch
            j, dJ = compFunc(shI, x)
            bDa[iF] = false
        end

        Jc[iF] = j
        dJc[iF] = mergeGradient(dJ)
    end

    jV = 0.0
    djV = zeros(size(dJc[1]))

    for iF in 1:nF
        wI = obj.weights[iF]
        jV += wI * Jc[iF]
        djV .+= wI * dJc[iF]
    end

    obj.isBatchDepleted = any(bDa)
    obj.shapeValues = Jc

    return jV, djV
end

function mergeGradient(dJ)
    if isa(dJ, Vector) && all(x -> hasproperty(x, :fValues), dJ)
        nDV = length(dJ)
        nDim1 = length(dJ[1].fValues)
        dJm = zeros(nDV * nDim1)

        for i in 1:nDV
            ind1 = 1 + nDim1 * (i - 1)
            ind2 = nDim1 + nDim1 * (i - 1)
            dJm[ind1:ind2] = dJ[i].fValues
        end

        return dJm
    elseif isa(dJ, AbstractVector{<:Number})
        return dJ
    else
        @warn "Unsupported input type. dJ should be a vector of structs or a numeric array."
        return zeros(1)
    end
end

end