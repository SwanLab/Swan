module Network

export Net, init_network, compute_output, forward_pass, backpropagation, compute_last_H, compute_gradient

using LinearAlgebra
using LearnableVariables

"""
    Net

Immutable struct describing the architecture and activation types of a neural network.
"""
struct Net
    neurons_per_layer::Vector{Int}
    n_layers::Int
    hu_type::String
    ou_type::String
    learnable_variables::LearnableVars
end

"""
    init_network(params)

Creates a `Net` object from parameter dictionary.
"""
function init_network(params::Dict{String, Any})
    hidden_layers = Vector{Int}(params["hiddenLayers"])
    data = params["data"]
    n_inputs = size(data.Xtrain, 2)
    n_outputs = data.nLabels

    neurons = [n_inputs; hidden_layers; n_outputs]
    n_layers = length(neurons)

    lv = init_learnable_vars(Dict(
        "neuronsPerLayer" => neurons,
        "nLayers" => n_layers
    ))

    return Net(neurons, n_layers, params["HUtype"], params["OUtype"], lv)
end

"""
    compute_output(net, X)

Computes network output for inputs X.
"""
function compute_output(net::Net, X::Matrix{Float64})
    a_values = forward_pass(net, X)
    return a_values[end]
end

"""
    forward_pass(net, X)

Performs a forward pass and returns (z_values, a_values), where:
- z_values[k] is pre-activation at layer k
- a_values[k] is activation at layer k
"""
function forward_pass(net::Net, X::Matrix{Float64})
    W, b = reshape_in_layer_form(net.learnable_variables)
    a_values = Vector{Matrix{Float64}}(undef, net.n_layers)
    a_values[1] = X

    for i in 2:net.n_layers
        z = a_values[i-1] * W[i-1] .+ b[i-1]'
        act_type = i == net.n_layers ? net.ou_type : net.hu_type
        a = activation(z, act_type)
        a_values[i] = a
    end

    return a_values
end

"""
    backpropagation(net, Y, dLF, a_values)

Computes the gradient of the loss wrt weights using backpropagation.
Returns a flattened gradient vector.
"""
function backpropagation(net::Net, Y::Matrix{Float64}, dLF::Matrix{Float64}, a_values::Vector{Matrix{Float64}})
    W, _ = reshape_in_layer_form(net.learnable_variables)
    npl = net.neurons_per_layer
    m = size(Y, 1)
    nL = net.n_layers

    deltag = Vector{Matrix{Float64}}(undef, nL)
    dcW = Vector{Matrix{Float64}}(undef, nL - 1)
    dcB = Vector{Vector{Float64}}(undef, nL - 1)

    for k in reverse(2:nL)
        z_k = a_values[k]
        act_type = k == nL ? net.ou_type : net.hu_type
        g_der = activation_derivative(z_k, act_type)

        if k == nL
            deltag[k] = dLF .* g_der
        else
            deltag[k] = deltag[k+1] * W[k]' .* g_der
        end

        dcW[k-1] = a_values[k-1]' * deltag[k] ./ m
        dcB[k-1] = vec(sum(deltag[k], dims=1)) ./ m
    end

    # Flatten dcW and dcB into gradient vector
    grad = Float64[]
    for i in 1:length(dcW)
        append!(grad, vec(dcW[i]))
        append!(grad, dcB[i])
    end
    return grad
end

"""
    compute_gradient(net, X)

Computes ∂a/∂x for final output layer wrt X.
"""
function compute_gradient(net::Net, X::Matrix{Float64})
    W, _ = reshape_in_layer_form(net.learnable_variables)
    a_values = forward_pass(net, X)
    nL = net.n_layers

    for k in reverse(1:nL-1)
        act_type = k+1 == nL ? net.ou_type : net.hu_type
        a_der = activation_derivative(a_values[k+1], act_type)
        A = Diagonal(vec(a_der))
        parDer = A * W[k]'
        grad = k+1 == nL ? parDer : grad * parDer
    end
    return grad
end

"""
    compute_last_H(net, X)

Returns the last hidden layer activation (before output layer).
"""
function compute_last_H(net::Net, X::Matrix{Float64})
    W, b = reshape_in_layer_form(net.learnable_variables)
    a = X
    for i in 1:(net.n_layers - 1)
        z = a * W[i] .+ b[i]'
        act_type = i == net.n_layers - 1 ? net.ou_type : net.hu_type
        a = activation(z, act_type)
    end
    return a
end

# --- Activation functions and their derivatives ---

function activation(z::Matrix{Float64}, type::String)
    if type == "sigmoid"
        return 1.0 ./ (1.0 .+ exp.(-z))
    elseif type == "tanh"
        return tanh.(z)
    elseif type == "ReLU"
        return max.(z, 0.0)
    elseif type == "linear"
        return z
    elseif type == "softmax"
        ez = exp.(z .- maximum(z, dims=2))
        return ez ./ sum(ez, dims=2)
    else
        error("Unknown activation type: $type")
    end
end

function activation_derivative(z::Matrix{Float64}, type::String)
    if type == "sigmoid"
        g = activation(z, type)
        return g .* (1 .- g)
    elseif type == "tanh"
        g = tanh.(z)
        return 1 .- g.^2
    elseif type == "ReLU"
        return Float64.(z .> 0.0)
    elseif type == "linear"
        return ones(size(z))
    elseif type == "softmax"
        g = activation(z, type)
        return g .* (1 .- g)  # Approximate
    else
        error("Unknown activation type: $type")
    end
end

end