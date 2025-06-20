module LearnableVariables

export LearnableVars, reshapeInLayerForm

using Random
using Distributions

struct LearnableVars
    neuronsPerLayer::Vector{Int}
    nLayers::Int
    thetavec::Vector{Float64}
end

# Constructor: builds and returns a LearnableVars object
function LearnableVars(s::Dict{String, Any})
    nPL = Vector{Int}(s["neuronsPerLayer"])
    nLayers = s["nLayers"]
    thetavec = computeInitialTheta(nPL, nLayers)
    return LearnableVars(nPL, nLayers, thetavec)
end

# Reshape theta into weights and biases for each layer
function reshapeInLayerForm(obj::LearnableVars)
    theta = obj.thetavec
    nPL = obj.neuronsPerLayer
    last = 1
    b = Vector{Vector{Float64}}(undef, obj.nLayers - 1)
    W = Vector{Matrix{Float64}}(undef, obj.nLayers - 1)

    for i in 2:obj.nLayers
        prevL = nPL[i-1]
        nextL = nPL[i]
        aux = prevL * nextL + nextL
        next = last + aux - 1
        thetaI = theta[last:next]
        b[i-1] = getB(thetaI, prevL, nextL)
        W[i-1] = getW(thetaI, prevL, nextL)
        last = next + 1
    end

    return W, b
end

function computeInitialTheta(nPL::Vector{Int}, nLayers::Int)
    th = Float64[]
    for i in 2:nLayers
        nextL = nPL[i]
        prevL = nPL[i-1]

        # Bias vector
        if i != nLayers
            getB = fill(0.1, nextL)
        else
            getB = fill(1 / nextL, nextL)
        end

        # Weight matrix flattened
        u = sqrt(6 / (prevL + nextL))
        getW = rand(Uniform(-u, u), prevL * nextL)

        append!(th, getW)
        append!(th, getB)
    end
    return th
end

function getW(thv::Vector{Float64}, prevL::Int, nextL::Int)
    aux = thv[1:prevL * nextL]
    return reshape(aux, (prevL, nextL))
end

function getB(thv::Vector{Float64}, prevL::Int, nextL::Int)
    return thv[prevL * nextL + 1:end]
end

end