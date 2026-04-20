module Network

export Net, NetworkBuffers, init_network, init_buffers,
       compute_output, forward_pass!, backpropagation!,
       compute_last_H, compute_gradient

using LinearAlgebra
using ..LearnableVariables

struct Net
    neurons_per_layer::Vector{Int}
    n_layers::Int
    hu_type::String
    ou_type::String
    learnable_variables::LearnableVars
end

"""
    NetworkBuffers

Buffers pré-alloués pour forward pass et backpropagation.
Créés une seule fois à l'init, réutilisés à chaque step.
tmp : buffer temporaire pour deltag[k+1] * W[k]' dans backprop.
"""
struct NetworkBuffers
    a_values :: Vector{Matrix{Float64}}
    deltag   :: Vector{Matrix{Float64}}
    dcW      :: Vector{Matrix{Float64}}
    dcB      :: Vector{Vector{Float64}}
    grad     :: Vector{Float64}
    tmp      :: Vector{Matrix{Float64}}   # NOUVEAU : évite l'allocation de tmp dans backprop
end

"""
    init_buffers(net, batch_size)

Alloue tous les buffers nécessaires pour un batch de taille donnée.
À appeler une seule fois depuis Adam.
"""
function init_buffers(net::Net, batch_size::Int)
    npl = net.neurons_per_layer
    nL  = net.n_layers

    a_values = Vector{Matrix{Float64}}(undef, nL)
    a_values[1] = zeros(batch_size, npl[1])
    for i in 2:nL
        a_values[i] = zeros(batch_size, npl[i])
    end

    deltag = Vector{Matrix{Float64}}(undef, nL)
    for k in 2:nL
        deltag[k] = zeros(batch_size, npl[k])
    end

    dcW = [zeros(npl[i], npl[i+1]) for i in 1:nL-1]
    dcB = [zeros(npl[i+1])          for i in 1:nL-1]
    grad = zeros(sum(npl[i]*npl[i+1] + npl[i+1] for i in 1:nL-1))

    # tmp[k] a la taille de deltag[k] : batch_size × npl[k]
    # utilisé pour stocker deltag[k+1] * W[k]' avant multiplication élément par élément
    tmp = Vector{Matrix{Float64}}(undef, nL)
    for k in 2:nL
        tmp[k] = zeros(batch_size, npl[k])
    end

    return NetworkBuffers(a_values, deltag, dcW, dcB, grad, tmp)
end

function init_network(params::Dict{String, Any})
    hidden_layers = Vector{Int}(params["hiddenLayers"])
    data          = params["data"]
    n_inputs      = size(data.Xtrain, 2)
    n_outputs     = data.n_labels

    neurons  = [n_inputs; hidden_layers; n_outputs]
    n_layers = length(neurons)

    lv = init_learnable_variables(Dict(
        "neuronsPerLayer" => neurons,
        "nLayers"         => n_layers
    ))

    return Net(neurons, n_layers, params["HUtype"], params["OUtype"], lv)
end

"""
    forward_pass!(buf, net, X, W, b)

Forward pass in-place. W et b sont passés depuis Adam — zéro reshape par step.
"""
function forward_pass!(buf::NetworkBuffers, net::Net,
                        X::AbstractMatrix{Float64},
                        W::Vector{<:AbstractMatrix{Float64}},
                        b::Vector{<:AbstractVector{Float64}})
    nL = net.n_layers
    buf.a_values[1] = X   # vue, pas de copie

    for i in 2:nL
        act_type = (i == nL) ? net.ou_type : net.hu_type
        mul!(buf.a_values[i], buf.a_values[i-1], W[i-1])
        buf.a_values[i] .+= b[i-1]'
        _activation!(buf.a_values[i], act_type)
    end
end

"""
    backpropagation!(buf, net, Y, dLF, W)

Backpropagation in-place. W passé depuis Adam — zéro reshape par step.
tmp pré-alloué dans buf évite l'allocation de la matrice intermédiaire.
"""
function backpropagation!(buf::NetworkBuffers, net::Net,
                           Y::AbstractMatrix{Float64},
                           dLF::Matrix{Float64},
                           W::Vector{<:AbstractMatrix{Float64}})
    npl = net.neurons_per_layer
    m   = size(Y, 1)
    nL  = net.n_layers

    for k in reverse(2:nL)
        act_type = (k == nL) ? net.ou_type : net.hu_type

        _activation_derivative!(buf.deltag[k], buf.a_values[k], act_type)

        if k == nL
            buf.deltag[k] .*= dLF
        else
            # AVANT : tmp = buf.deltag[k+1] * W[k]'  → allocation
            # APRÈS : mul! in-place dans buf.tmp[k]   → zéro allocation
            mul!(buf.tmp[k], buf.deltag[k+1], W[k]')
            buf.deltag[k] .*= buf.tmp[k]
        end

        mul!(buf.dcW[k-1], buf.a_values[k-1]', buf.deltag[k])
        buf.dcW[k-1] ./= m

        sum!(reshape(buf.dcB[k-1], 1, :), buf.deltag[k])
        buf.dcB[k-1] ./= m
    end

    # Aplatissement in-place dans buf.grad
    offset = 1
    for i in 1:nL-1
        nW = npl[i] * npl[i+1]
        buf.grad[offset : offset + nW - 1]                  .= vec(buf.dcW[i])
        buf.grad[offset + nW : offset + nW + npl[i+1] - 1] .= buf.dcB[i]
        offset += nW + npl[i+1]
    end
end

# ---------------------------------------------------------------------------
# Fonctions non-buffered — utilisées uniquement hors boucle d'entraînement
# ---------------------------------------------------------------------------

function compute_output(net::Net, X::Matrix{Float64}, θ::Vector{Float64})
    W, b = reshape_in_layer_form(update_thetavec(net.learnable_variables, θ))
    nL   = net.n_layers
    a    = X
    for i in 2:nL
        act_type = (i == nL) ? net.ou_type : net.hu_type
        a = _activation_alloc(a * W[i-1] .+ b[i-1]', act_type)
    end
    return a
end

function compute_last_H(net::Net, X::Matrix{Float64}, θ::Vector{Float64})
    W, b = reshape_in_layer_form(update_thetavec(net.learnable_variables, θ))
    a    = X
    for i in 1:(net.n_layers - 2)
        a = _activation_alloc(a * W[i] .+ b[i]', net.hu_type)
    end
    return a
end

function compute_gradient(net::Net, X::Matrix{Float64}, θ::Vector{Float64})
    W, _     = reshape_in_layer_form(update_thetavec(net.learnable_variables, θ))
    nL       = net.n_layers
    a_values = Vector{Matrix{Float64}}(undef, nL)
    a_values[1] = X
    for i in 2:nL
        act_type    = (i == nL) ? net.ou_type : net.hu_type
        a_values[i] = _activation_alloc(a_values[i-1] * W[i-1] .+ b[i-1]', act_type)
    end
    grad = Diagonal(vec(_activation_deriv_alloc(a_values[nL], net.ou_type))) * W[nL-1]'
    for k in reverse(1:nL-2)
        g_der = _activation_deriv_alloc(a_values[k+1], net.hu_type)
        grad  = grad * (Diagonal(vec(g_der)) * W[k]')
    end
    return grad
end

# ---------------------------------------------------------------------------
# Activations in-place
# ---------------------------------------------------------------------------

function _activation!(a::Matrix{Float64}, type::String)
    if type == "ReLU"
        @. a = max(a, 0.0)
    elseif type == "linear"
        # rien à faire
    elseif type == "sigmoid"
        @. a = 1.0 / (1.0 + exp(-a))
    elseif type == "tanh"
        @. a = tanh(a)
    else
        error("Type d'activation inconnu : $type")
    end
end

function _activation_derivative!(out::Matrix{Float64}, a::Matrix{Float64}, type::String)
    if type == "ReLU"
        @. out = Float64(a > 0.0)
    elseif type == "linear"
        fill!(out, 1.0)
    elseif type == "sigmoid"
        @. out = a * (1.0 - a)
    elseif type == "tanh"
        @. out = 1.0 - tanh(a)^2
    else
        error("Type d'activation inconnu : $type")
    end
end

# Versions allouantes — uniquement pour compute_output / compute_gradient
function _activation_alloc(z::Matrix{Float64}, type::String)
    if type == "ReLU";       return max.(z, 0.0)
    elseif type == "linear"; return z
    elseif type == "sigmoid"; return @. 1.0 / (1.0 + exp(-z))
    elseif type == "tanh";   return tanh.(z)
    else error("Type d'activation inconnu : $type"); end
end

function _activation_deriv_alloc(a::Matrix{Float64}, type::String)
    if type == "ReLU";       return Float64.(a .> 0.0)
    elseif type == "linear"; return ones(size(a))
    elseif type == "sigmoid"; return @. a * (1.0 - a)
    elseif type == "tanh";   return @. 1.0 - tanh(a)^2
    else error("Type d'activation inconnu : $type"); end
end

end