module Network

export Net, computeYOut, backprop, networkGradient, computeLastH, getLearnableVariables

include("../LearnableVariables/LearnableVariables.jl")
using .LearnableVariables 

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
    aValues::Union{Vector{Any}, Nothing}
    deltag::Union{Vector{Any}, Nothing}
end

function Net(cParams::Dict{String, Any})
    data = cParams["data"]
    hiddenLayers = Vector{Int}(cParams["hiddenLayers"])
    nFeatures = data["nFeatures"]

    # Necessary for Matlab communication only, conditional for robustness
    Xraw = data["Xtrain"]
    if ndims(Xraw) == 2
        # Already a matrix, keep as is
        Xtrain = Xraw
    else
        # Assume it's a vector of vectors, concatenate and transpose
        Xtrain = hcat(Xraw...)'
    end
    data["Xtrain"] = Xtrain


    nPolyFeatures = size(data["Xtrain"], 2)
    nLabels = data["nLabels"]
    HUtype = cParams["HUtype"]
    OUtype = cParams["OUtype"]
    neuronsPerLayer = [nPolyFeatures; hiddenLayers; nLabels]
    nLayers = length(neuronsPerLayer)

    learnableParams = Dict(
        "neuronsPerLayer" => neuronsPerLayer,
        "nLayers" => nLayers
    )

    # Pass thetavec if available in full args (maybe only Matlab necessary)
    if haskey(cParams, "thetavec")
        println("Net received thetavec, passing it to LearnableVars")
        learnableParams["thetavec"] = cParams["thetavec"]
    end

    learnableVariables = LearnableVars(learnableParams) # equivalent to createLearnableVariables

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
        nothing, # aValues
        nothing  # deltag
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
    m = size(Yb, 1)

    #obj.deltag = fill(nothing, nLy)
    dcW = Vector{Any}(undef, nLy - 1)
    dcB = Vector{Any}(undef, nLy - 1)
    obj.deltag = Vector{Any}(undef, nLy)
    obj.deltag[1] = zeros(0)  # Match MATLAB’s {0x0 double}

    for k in reverse(2:nLy)
        _, g_der = actFCN(obj, a[k], k)
        if k == nLy
            obj.deltag[k] = dLF .* g_der
        else
            #obj.deltag[k] = (W[k] * obj.deltag[k+1]')' .* g_der
            obj.deltag[k] = (obj.deltag[k+1] * W[k]') .* g_der
        end
        dcW[k-1] = (1 / m) * (a[k-1]' * obj.deltag[k])
        dcB[k-1] = vec(sum(obj.deltag[k], dims=1)) / m
        #dcB[k-1] = (1 / m) * sum(obj.deltag[k], dims=1)
    end
    #println(obj.deltag)
    dc = Float64[]
    for i in 2:nLy
        aux1 = vcat(vec(dcW[i-1]), vec(dcB[i-1]))
        append!(dc, aux1)
    end
    return dc
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
    a = Vector{Any}(undef, nLy)
    a[1] = X
    for i in 2:nLy
        g_prev = a[i-1]
        Wi = W[i-1]
        bi = b[i-1]
        h = hypothesisfunction(g_prev, Wi, bi)
        g, _ = actFCN(obj, h, i)
        a[i] = g
    end
    obj.aValues = a
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
