module Optimizer

export AdamStruct, SGDStruct, init_Adam, init_SGD, compute, plot_cost_func

using ..CostNN
using ..Network
using ..LearnableVariables
using LinearAlgebra
using Random
using Plots
using Printf

# ---------------------------------------------------------------------------
# Structs
# ---------------------------------------------------------------------------

struct AdamStruct
    cost_function :: CostNN.CostNNStruct
    α             :: Float64
    β1            :: Float64
    β2            :: Float64
    ε             :: Float64
    batch_size    :: Int
    max_epochs    :: Int
    m             :: Vector{Float64}
    v             :: Vector{Float64}
    t             :: Int
    fplot         :: Vector{Float64}
end

struct SGDStruct
    cost_function :: CostNN.CostNNStruct
    α             :: Float64
    batch_size    :: Int
    max_epochs    :: Int
    fplot         :: Vector{Float64}
end

# ---------------------------------------------------------------------------
# Initialiseurs
# ---------------------------------------------------------------------------

function init_Adam(params::Dict{String, Any})
    cost_func  = params["costFunc"]
    lv         = params["designVariable"]
    n_params   = length(lv.thetavec)
    n_samples  = size(cost_func.shapeFunctions[1].data.Xtrain, 1)
    batch_size = min(n_samples, 200)

    return AdamStruct(
        cost_func,
        Float64(params["learningRate"]),
        Float64(get(params, "β1",  0.9)),
        Float64(get(params, "β2",  0.999)),
        Float64(get(params, "ε",   1e-8)),
        batch_size,
        params["maxEpochs"],
        zeros(n_params),
        zeros(n_params),
        0,
        zeros(params["maxEpochs"])
    )
end

function init_SGD(params::Dict{String, Any})
    cost_func  = params["costFunc"]
    n_samples  = size(cost_func.shapeFunctions[1].data.Xtrain, 1)
    batch_size = min(n_samples, 200)

    return SGDStruct(
        cost_func,
        Float64(params["learningRate"]),
        batch_size,
        params["maxEpochs"],
        zeros(params["maxEpochs"])
    )
end

# ---------------------------------------------------------------------------
# Boucle commune — extraite pour éviter la duplication
# ---------------------------------------------------------------------------

"""
    _training_loop(cost_function, buf, θ, fplot, batch_size, max_epochs, update_fn!)

Boucle d'entraînement partagée entre Adam et SGD.
`update_fn!` est la fonction de mise à jour spécifique à l'optimiseur :
  - reçoit (θ, grad, state) et retourne state mis à jour
  - c'est le seul endroit qui diffère entre Adam et SGD
"""
function _training_loop(cost_function, net, Xtrain, Ytrain,
                         buf, θ, fplot, batch_size, max_epochs,
                         w2, update_fn!, state, start_time)

    n_samples = size(Xtrain, 1)
    n_batches = max(1, fld(n_samples, batch_size))
    n_outputs = size(Ytrain, 2)
    n_inputs  = size(Xtrain, 2)

    dLF    = zeros(batch_size, n_outputs)
    order  = Vector{Int}(undef, n_samples)
    Xb_buf = Matrix{Float64}(undef, batch_size, n_inputs)
    Yb_buf = Matrix{Float64}(undef, batch_size, n_outputs)

    for epoch in 1:max_epochs
        randperm!(order)

        lv   = update_thetavec(net.learnable_variables, θ)
        W, b = reshape_in_layer_form(lv)

        epoch_cost = 0.0

        for batch in 1:n_batches
            i_start   = (batch - 1) * batch_size + 1
            i_end     = min(batch * batch_size, n_samples)
            actual_bs = i_end - i_start + 1

            # Copie explicite dans buffers contigus
            @inbounds for j in 1:n_inputs
                for (bi, si) in enumerate(@view order[i_start:i_end])
                    Xb_buf[bi, j] = Xtrain[si, j]
                end
            end
            @inbounds for j in 1:n_outputs
                for (bi, si) in enumerate(@view order[i_start:i_end])
                    Yb_buf[bi, j] = Ytrain[si, j]
                end
            end

            Xb = @view Xb_buf[1:actual_bs, :]
            Yb = @view Yb_buf[1:actual_bs, :]

            Network.forward_pass!(buf, net, Xb, W, b)

            J   = _compute_cost_and_dLF!(dLF, cost_function, Yb, buf.a_values[end])
            J  += w2 * 0.5 * dot(θ, θ)

            Network.backpropagation!(buf, net, Yb, dLF, W)
            @. buf.grad += w2 * θ

            # Mise à jour spécifique à l'optimiseur — seul point de divergence
            state = update_fn!(θ, buf.grad, state)

            epoch_cost += J
        end

        # W et b périmés après update θ — recalculés à la prochaine époque
        fplot[epoch] = epoch_cost / n_batches

        if epoch % 100 == 0 || epoch == 1
            @printf("Époque %4d / %d  |  coût = %.6g  |  %.1f s\n",
                    epoch, max_epochs, fplot[epoch],
                    round(time() - start_time, digits=1))
        end
    end

    return state
end

# ---------------------------------------------------------------------------
# compute — dispatch multiple : Julia choisit selon le type de l'optimiseur
# ---------------------------------------------------------------------------

function compute(adam::AdamStruct, θ::Vector{Float64})
    start_time = time()

    lf     = adam.cost_function.shapeFunctions[1]
    net    = lf.network
    buf    = Network.init_buffers(net, adam.batch_size)
    fplot  = copy(adam.fplot)
    w2     = adam.cost_function.weights[2]

    # État Adam : (m, v, t, hyperparamètres)
    α, β1, β2, ε = adam.α, adam.β1, adam.β2, adam.ε
    state = (m=copy(adam.m), v=copy(adam.v), t=adam.t,
             α=α, β1=β1, β2=β2, ε=ε)

    # Fonction de mise à jour Adam — capture les hyperparamètres par fermeture
    function adam_update!(θ, grad, s)
        t_new = s.t + 1
        @. s.m = s.β1 * s.m + (1.0 - s.β1) * grad
        @. s.v = s.β2 * s.v + (1.0 - s.β2) * grad^2
        α_t = s.α * sqrt(1.0 - s.β2^t_new) / (1.0 - s.β1^t_new)
        @. θ -= α_t * s.m / (sqrt(s.v) + s.ε)
        return (m=s.m, v=s.v, t=t_new, α=s.α, β1=s.β1, β2=s.β2, ε=s.ε)
    end

    state = _training_loop(
        adam.cost_function, net, lf.data.Xtrain, lf.data.Ytrain,
        buf, θ, fplot, adam.batch_size, adam.max_epochs,
        w2, adam_update!, state, start_time
    )

    println("Entraînement terminé. Durée totale : $(round(time()-start_time, digits=2)) s")

    adam_final = AdamStruct(adam.cost_function, α, β1, β2, ε,
                             adam.batch_size, adam.max_epochs,
                             state.m, state.v, state.t, fplot)
    return adam_final, θ
end

function compute(sgd::SGDStruct, θ::Vector{Float64})
    start_time = time()

    lf    = sgd.cost_function.shapeFunctions[1]
    net   = lf.network
    buf   = Network.init_buffers(net, sgd.batch_size)
    fplot = copy(sgd.fplot)
    w2    = sgd.cost_function.weights[2]
    α     = sgd.α

    # État SGD : juste le learning rate
    state = (α=α,)

    # Fonction de mise à jour SGD — descente de gradient simple
    function sgd_update!(θ, grad, s)
        @. θ -= s.α * grad
        return s   # état inchangé pour SGD statique
    end

    state = _training_loop(
        sgd.cost_function, net, lf.data.Xtrain, lf.data.Ytrain,
        buf, θ, fplot, sgd.batch_size, sgd.max_epochs,
        w2, sgd_update!, state, start_time
    )

    println("Entraînement terminé. Durée totale : $(round(time()-start_time, digits=2)) s")

    sgd_final = SGDStruct(sgd.cost_function, α, sgd.batch_size, sgd.max_epochs, fplot)
    return sgd_final, θ
end

# ---------------------------------------------------------------------------
# Utilitaires communs
# ---------------------------------------------------------------------------

function _compute_cost_and_dLF!(dLF::Matrix{Float64},
                                  cost::CostNN.CostNNStruct,
                                  Y::AbstractMatrix{Float64},
                                  Ŷ::AbstractMatrix{Float64})
    w1   = cost.weights[1]
    diff = Ŷ .- Y
    J    = w1 * sqrt(sum(diff .^ 2))
    @. dLF = w1 * diff
    return J
end

function plot_cost_func(opt::AdamStruct)
    _plot_fplot(opt.fplot, "Convergence Adam")
end

function plot_cost_func(opt::SGDStruct)
    _plot_fplot(opt.fplot, "Convergence SGD")
end

function _plot_fplot(fplot::Vector{Float64}, title_str::String)
    valid = fplot[fplot .> 0.0]
    isempty(valid) && (println("Aucune donnée de coût."); return)
    plt = plot(1:length(valid), valid;
               xlabel="Époques", ylabel="Coût", title=title_str,
               linewidth=1.8, legend=false, grid=true, yscale=:log10)
    display(plt)
    return plt
end

end