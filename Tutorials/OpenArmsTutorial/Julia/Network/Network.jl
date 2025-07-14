module Network

export Net, computeYOut, backprop, networkGradient, computeLastH, getLearnableVariables

include("../LearnableVariables/LearnableVariables.jl")
using .LearnableVariables
using LinearAlgebra

mutable struct Net
    hiddenLayers::Vector{Int}
    nFeatures::Int
    nPolyFeatures::Int
    nLabels::Int
    neuronsPerLayer::Vector{Int}
    nLayers::Int
    HUtype::String
    OUtype::String
    learnableVariables::LearnableVars
    aValues::Vector{Matrix{Float64}}
    deltag::Vector{Matrix{Float64}}
    dcW::Vector{Matrix{Float64}}
    dcB::Vector{Vector{Float64}}
    gradient::Vector{Float64}
end

function Net(cParams::Dict{String, Any})
    data = cParams["data"]
    hiddenLayers = Vector{Int}(cParams["hiddenLayers"])
    nFeatures = data.nFeatures

    nPolyFeatures = size(data.Xtrain, 2)
    nLabels = data.nLabels
    HUtype = cParams["HUtype"]
    OUtype = cParams["OUtype"]
    neuronsPerLayer = [nPolyFeatures; hiddenLayers; nLabels]
    nLayers = length(neuronsPerLayer)

    learnableParams = Dict(
        "neuronsPerLayer" => neuronsPerLayer,
        "nLayers" => nLayers
    )

    learnableVariables = LearnableVars(learnableParams) # equivalent to createLearnableVariables

    aValues = [zeros(0,0) for _ in 1:nLayers]
    
    # Initialization of the backprop parameters
    deltag = [zeros(size(aValues[i])) for i in 1:nLayers]
    dcW = [zeros(neuronsPerLayer[k-1], neuronsPerLayer[k]) for k in 2:nLayers]
    dcB = [zeros(neuronsPerLayer[i+1]) for i in 1:nLayers-1]
    gradSize = sum(length.(dcW)) + sum(length.(dcB)) # Reshape whole gradient in one vector
    gradient = zeros(gradSize)

    return Net(
        hiddenLayers,
        nFeatures,
        nPolyFeatures,
        nLabels,
        neuronsPerLayer,
        nLayers,
        HUtype,
        OUtype,
        learnableVariables,
        aValues, # aValues
        deltag,  # deltag
        dcW,
        dcB,
        gradient
    )
end

function computeYOut(obj::Net, Xb::Matrix{Float64})
    computeAvalues(obj, Xb)
    return obj.aValues[end]
end

function backprop(obj::Net, Yb::Matrix{Float64}, dLF::Matrix{Float64})
    W, _ = reshapeInLayerForm(obj.learnableVariables)
    a = obj.aValues
    nPl = obj.neuronsPerLayer
    nLy = obj.nLayers
    m = size(Yb, 1) # Batch size

    #obj.deltag = fill(nothing, nLy)
    #=
    dcW = Vector{Any}(undef, nLy - 1)
    dcB = Vector{Any}(undef, nLy - 1)
    =#
    #dcW = [zeros(nPl[k-1], nPl[k]) for k in 2:nLy]
    #dcB = [zeros(nPl[i+1]) for i in 1:nLy-1]
    #obj.deltag = [zeros(m, nPl[i]) for i in 1:nLy]

    # Ensure deltag has correct shape and size
    for i in 1:nLy
        if size(obj.deltag[i], 1) != m || size(obj.deltag[i], 2) != nPl[i]
            obj.deltag[i] = zeros(m, nPl[i])
        end
    end

    for k in reverse(2:nLy)
        _, g_der = actFCN(obj, a[k], k)
        if k == nLy
            obj.deltag[k] .= dLF .* g_der
        else
            mul!(obj.deltag[k], obj.deltag[k+1], W[k]')
            obj.deltag[k] .*= g_der
        end

        mul!(obj.dcW[k-1], a[k-1]', obj.deltag[k])
        obj.dcW[k-1] ./= m
        obj.dcB[k-1] .= sum(obj.deltag[k], dims=1)' ./ m
    end
    
    # Different to Matlab's code for performance
    # Start filling the gradient vector from index 1
    offset = 0
    # Loop through each layer's gradients
    for i in 1:length(obj.dcW)
        dw = obj.dcW[i]  # weight gradient matrix of layer i
        db = obj.dcB[i]  # bias gradient vector of layer i

        # Flatten weight matrix and copy it into the correct position in gradient
        copyto!(obj.gradient, offset + 1, vec(dw), 1, length(vec(dw)))
        offset += length(vec(dw))  # Update offset for next copy

        # Copy bias vector into the gradient vector
        copyto!(obj.gradient, offset + 1, db, 1, length(db))
        offset += length(db)  # Update offset again
    end
    return obj.gradient
    
end

function networkGradient(obj::Net, X::Matrix{Float64})
    computeAvalues(obj, X)
    W, _ = reshapeInLayerForm(obj.learnableVariables)
    a = obj.aValues
    nLy = obj.nLayers

    for k in reverse(1:(nLy-1))
        _, a_der = actFCN(obj, a[k+1], k+1)
        a_der_diag = Diagonal(vec(a_der))
        parDer = a_der_diag * W[k]'
        if k == nLy - 1
            grad = parDer
        else
            grad = grad * parDer
        end
    end
    return grad
end

function computeLastH(obj::Net, X::Matrix{Float64})
    nLy = obj.nLayers
    W, b = reshapeInLayerForm(obj.learnableVariables)

    h = hypothesisfunction(X, W[1], b[1])
    g, _ = actFCN(obj, h, 2)

    for i in 2:(nLy - 1)
        h = hypothesisfunction(g, W[i], b[i])
        g, _ = actFCN(obj, h, i + 1)
    end

    return g
end

function getLearnableVariables(obj::Net)
    return obj.learnableVariables
end

function computeAvalues(obj::Net, X::Matrix{Float64})
    W, b = reshapeInLayerForm(obj.learnableVariables)
    nLy = obj.nLayers
    obj.aValues[1] = X
    for i in 2:nLy
        g_prev = obj.aValues[i-1]
        Wi = W[i-1]
        bi = b[i-1]
        h = hypothesisfunction(g_prev, Wi, bi)
        g, _ = actFCN(obj, h, i)
        obj.aValues[i] = g
    end
end

function actFCN(obj::Net, z::Matrix{Float64}, k::Int)
    nLy = obj.nLayers
    type = k == nLy ? obj.OUtype : obj.HUtype
    if type == "sigmoid"
        g = 1 ./ (1 .+ exp.(-z))
        g_der = g .* (1 .- g)
    elseif type == "ReLU"
        g = max.(z, 0.0)
        g_der = Float64.(z .> 0.0)
    elseif type == "tanh"
        g = tanh.(z)
        g_der = 1 .- g.^2
    elseif type == "softmax"
        ez = exp.(z .- maximum(z, dims=2))  # for numerical stability
        g = ez ./ sum(ez, dims=2)
        g_der = g .* (1 .- g)  # note: true derivative is more complex
    elseif type == "linear"
        g = z
        g_der = ones(size(z))
    else
        error("$type is not a valid activation function")
    end
    return g, g_der
end

function hypothesisfunction(X::Matrix{Float64}, W::Matrix{Float64}, b::Vector{Float64})
    return X * W .+ b'
end

end
