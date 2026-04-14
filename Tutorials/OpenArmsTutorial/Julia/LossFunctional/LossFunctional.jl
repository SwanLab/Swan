module LossFunctional

export LossFunctionalStruct, computeFunctionAndGradient, computeStochasticCostAndGradient, getTestError, computeBatchSize


using ..Network
using ..LearnableVariables
using ..Data
using Random
using Distributions

mutable struct LossFunctionalStruct
    iBatch::Int
    order::Vector{Int}
    nBatches::Int

    costType::String
    designVariable::LearnableVars   
    network::Any
    data::DataStruct
end

function LossFunctionalStruct(cParams::Dict{String, Any})
    obj = LossFunctionalStruct(
        1,                              # iBatch starts at 1
        Int[],                         # order
        0,                             # nBatches
        cParams["costType"],
        cParams["designVariable"],
        cParams["network"],
        cParams["data"]
    )
    computeNumberOfBatchesAndOrder!(obj)
    return obj
end

function computeFunctionAndGradient(obj::LossFunctionalStruct, x::Vector{Float64})
    obj.designVariable.thetavec = x
    Xb = obj.data.Xtrain
    Yb = obj.data.Ytrain
    yOut = Network.computeYOut(obj.network, Xb)
    j  = computeCost(obj, yOut, Yb)
    dj = computeGradient(obj, yOut, Yb)
    return j, dj
end

function computeStochasticCostAndGradient(obj::LossFunctionalStruct, x::Vector{Float64}, moveBatch::Bool)
    obj.designVariable.thetavec = x
    Xt = obj.data.Xtrain
    Yt = obj.data.Ytrain
    Xb, Yb = updateSampledDataSet(obj, Xt, Yt, obj.iBatch)
    yOut = Network.computeYOut(obj.network, Xb)
    j = computeCost(obj, yOut, Yb)
    dj = computeGradient(obj, yOut, Yb)
    obj.iBatch = updateBatchCounter(obj, obj.iBatch, moveBatch)
    isBD = isBatchDepleted(obj, obj.iBatch, moveBatch)
    return j, dj, isBD
end


function getTestError(obj::LossFunctionalStruct)
    Xtest = obj.data.Xtest
    Ytest = obj.data.Ytest

    H = Network.computeLastH(obj.network, Xtest)  # Should return matrix of predictions
    Ypred = map(row -> findmax(row)[2], eachrow(H))      # Row-wise argmax (like max(..., [], 2))

    Ytarget = map(row -> findmax(row)[2], eachrow(Ytest))

    mismatches = sum(Ypred .!= Ytarget)             #number of incorrect prediction
    testError = mismatches / length(Ytarget)        # % of wrong prediction

    return testError
end

function computeCost(obj::LossFunctionalStruct, yOut, Yb)       #get the value of the cost
    J, _ = LossFunction(obj, Yb, yOut)
    return J
end

function computeGradient(obj::LossFunctionalStruct, yOut, Yb)   #get the value of the gradient for all the weights and bias
    _, dLF = LossFunction(obj, Yb, yOut)
    dj = Network.backprop(obj.network, Yb, dLF)
    return dj
end

function LossFunction(obj::LossFunctionalStruct, y::Matrix{Float64}, yOut::Matrix{Float64})   #compute the cost function and gradient between Ypred and Yreel
    type = obj.costType
    yp = yOut .- 1e-11          #avoid dividing by zero

    if type == "-loglikelihood"                     # for classification
        c = sum((1 .- y) .* (-log.(1 .- yp)) .+ y .* (-log.(yp)), dims=2)
        J = mean(c)
        gc = (yp .- y) ./ (yp .* (1 .- yp))

    elseif type == "L2"                             # for regression
        c = (yp .- y).^2
        J = sqrt(sum(c))
        gc = yp .- y

    else
        error("Loss type '$type' is not valid")
    end

    return J, gc
end

function computeNumberOfBatchesAndOrder!(obj::LossFunctionalStruct)     #gives the order of the samples and the nbr of batch to LossFunctionalStruct
    Xtrain = obj.data.Xtrain
    nD = size(Xtrain, 1)                    #number of samples
    batchSize = computeBatchSize(obj)
    obj.nBatches = fld(nD, batchSize)      #number of batch that we can do with all the data

    if obj.nBatches ≤ 1                    #mix the samples if nbr of batch >1 and if number of batch = 1 there is only one big batch without any mix
        obj.order = collect(1:nD)
        obj.nBatches = 1
    else
        obj.order = randperm(nD)
    end
end

function isBatchDepleted(obj::LossFunctionalStruct, iBatch::Int, moveBatch::Bool)       #says if you passed through all the batch
    return iBatch == obj.nBatches && moveBatch
end

function updateBatchCounter(obj::LossFunctionalStruct, iB::Int, moveBatch::Bool)        #says which batch we have to use
    if iB < obj.nBatches && moveBatch
        return iB + 1
    elseif iB == obj.nBatches && moveBatch          #starts an new epoch
        return 1
    else
        return iB
    end
end

function computeBatchSize(obj::LossFunctionalStruct)        # batch size = or lower than 200 samples
    n = size(obj.data.Xtrain, 1)
    return n > 200 ? 200 : n
end

function updateSampledDataSet(obj::LossFunctionalStruct, Xl::Matrix{Float64}, Yl::Matrix{Float64}, iBatch::Int)    #gives the index of the mixed samples contain in a batch
    batchSize = computeBatchSize(obj)
    startIdx = (iBatch - 1) * batchSize + 1
    endIdx = min(iBatch * batchSize, size(Xl, 1)) 
    idx = obj.order[startIdx:endIdx]
    return Xl[idx, :], Yl[idx, :]
end

end
