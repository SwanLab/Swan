module LearnableVariables

export LearnableVars, init_learnable_variables, reshape_in_layer_form, set_theta, update_thetavec

using Random
using Distributions

"""
    LearnableVars

Immutable struct to hold the network architecture and parameters in vectorized form.
"""
struct LearnableVars
    neurons_per_layer::Vector{Int}
    n_layers::Int
    thetavec::Vector{Float64}
end

"""
    init_learnable_variables(params)

Create a `LearnableVars` object from a layer structure and optionally a `thetavec`.
"""
function init_learnable_variables(params::Dict{String, Any})
    npl = Vector{Int}(params["neuronsPerLayer"])
    n_layers = params["nLayers"]
    thetavec = haskey(params, "thetavec") ?
        Float64.(params["thetavec"]) :
        _compute_initial_theta(npl)
    return LearnableVars(npl, n_layers, thetavec)
end

"""
    reshape_in_layer_form(lv)

Given a `LearnableVars` object, returns (W, b) as a vector of matrices and vectors.
"""
function reshape_in_layer_form(lv::LearnableVars)
    θ = lv.thetavec
    npl = lv.neurons_per_layer
    W = Vector{Matrix{Float64}}(undef, lv.n_layers - 1)
    b = Vector{Vector{Float64}}(undef, lv.n_layers - 1)

    offset = 1
    for i in 2:lv.n_layers
        in_size, out_size = npl[i-1], npl[i]
        total = in_size * out_size + out_size
        slice = θ[offset : offset + total - 1]
        W[i-1] = reshape(slice[1:in_size*out_size], in_size, out_size)
        b[i-1] = slice[in_size*out_size+1:end]
        offset += total
    end
    return W, b
end

function update_thetavec(lv::LearnableVars, θ_new::Vector{Float64})
    return LearnableVars(lv.neurons_per_layer, lv.n_layers, θ_new)
end

"""
    compute_initial_theta(npl)

Xavier initialization of thetavec from layer sizes.
"""
function _compute_initial_theta(npl::Vector{Int})
    θ = Float64[]
    for i in 2:length(npl)
        in_dim, out_dim = npl[i-1], npl[i]
        u = sqrt(6 / (in_dim + out_dim))
        W = rand(Uniform(-u, u), in_dim * out_dim)
        b = i != length(npl) ? fill(0.1, out_dim) : fill(1 / out_dim, out_dim)
        append!(θ, W)
        append!(θ, b)
    end
    return θ
end

function set_theta(lv::LearnableVars, θ::Vector{Float64})
    return LearnableVars(lv.neurons_per_layer, lv.n_layers, θ)
end

end
