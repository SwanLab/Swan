module LearnableVariables

export LearnableVars, init_learnable_variables, reshape_in_layer_form, update_thetavec

using Random
using Distributions

"""
    LearnableVars

Immutable struct containing the network architecture and parameters in vectorized form.
"""
struct LearnableVars
    neurons_per_layer::Vector{Int}
    n_layers::Int
    thetavec::Vector{Float64}
end

"""
    init_learnable_variables(params)

Creates a `LearnableVars` object from a parameter dictionary.
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

Returns (W, b) as vectors of matrices/vectors, without memory copying (using @view).
"""
function reshape_in_layer_form(lv::LearnableVars)
    θ   = lv.thetavec
    npl = lv.neurons_per_layer
    n   = lv.n_layers - 1
    W   = Vector{Matrix{Float64}}(undef, n)
    b   = Vector{Vector{Float64}}(undef, n)

    offset = 1
    for i in 1:n
        in_size  = npl[i]
        out_size = npl[i+1]
        nW = in_size * out_size

        W[i] = reshape(@view(θ[offset : offset + nW - 1]), in_size, out_size)
        b[i] = @view(θ[offset + nW : offset + nW + out_size - 1])

        offset += nW + out_size
    end
    return W, b
end

"""
    update_thetavec(lv, θ_new)

Returns a new `LearnableVars` instance with `θ_new` as parameters.
"""
function update_thetavec(lv::LearnableVars, θ_new::Vector{Float64})
    return LearnableVars(lv.neurons_per_layer, lv.n_layers, θ_new)
end

# --- Initialization ---

function _compute_initial_theta(npl::Vector{Int})
    # Pre-compute total size to avoid reallocations
    total = sum(npl[i-1] * npl[i] + npl[i] for i in 2:length(npl))
    θ = Vector{Float64}(undef, total)

    offset = 1
    for i in 2:length(npl)
        in_dim, out_dim = npl[i-1], npl[i]
        nW = in_dim * out_dim

        # Xavier uniform initialization
        u = sqrt(6.0 / (in_dim + out_dim))
        θ[offset : offset + nW - 1] .= rand(Uniform(-u, u), nW)

        # Bias initialization: 0.1 for hidden layers, 1/out_dim for output layer
        b_val = (i != length(npl)) ? 0.1 : 1.0 / out_dim
        θ[offset + nW : offset + nW + out_dim - 1] .= b_val

        offset += nW + out_dim
    end
    return θ
end

end