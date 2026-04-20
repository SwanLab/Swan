module Adam

export AdamStruct, init_Adam, compute, plot_cost_func

using ..CostNN
using ..Network
using LinearAlgebra
using Random
using Plots
using Printf

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

function init_Adam(params::Dict{String, Any})
    cost_func  = params["costFunc"]
    lv         = params["designVariable"]
    n_params   = length(lv.thetavec)
    n_samples  = size(cost_func.shapeFunctions[1].data.Xtrain, 1)
    batch_size = min(n_samples, 200)

    return AdamStruct(
        cost_func,
        Float64(params["learningRate"]),
        Float64(get(params, "β1", 0.9)),
        Float64(get(params, "β2", 0.999)),
        Float64(get(params, "ε",  1e-8)),
        batch_size,
        params["maxEpochs"],
        zeros(n_params),
        zeros(n_params),
        0,
        zeros(params["maxEpochs"])
    )
end

function compute(adam::AdamStruct, θ::Vector{Float64})
    start_time = time()

    lf         = adam.cost_function.shapeFunctions[1]
    Xtrain     = lf.data.Xtrain
    Ytrain     = lf.data.Ytrain
    net        = lf.network
    n_samples  = size(Xtrain, 1)
    batch_size = adam.batch_size
    n_batches  = max(1, fld(n_samples, batch_size))

    # Pré-allocations — toutes hors boucle
    buf   = Network.init_buffers(net, batch_size)
    m     = copy(adam.m)
    v     = copy(adam.v)
    t     = adam.t
    fplot = copy(adam.fplot)
    α, β1, β2, ε = adam.α, adam.β1, adam.β2, adam.ε
    w2    = adam.cost_function.weights[2]
    dLF   = zeros(batch_size, size(Ytrain, 2))
    order = Vector{Int}(undef, n_samples)   # NOUVEAU : pré-alloué, rempli par randperm!

    for epoch in 1:adam.max_epochs

        # AVANT : order = randperm(n_samples)  → allocation à chaque époque
        # APRÈS : randperm!(order)             → in-place, zéro allocation
        randperm!(order)

        # W, b extraits une seule fois par époque depuis θ courant
        # AVANT : reshape_in_layer_form appelé dans forward_pass! ET backpropagation!
        # APRÈS : un seul appel par époque, W et b partagés
        lv   = update_thetavec(net.learnable_variables, θ)
        W, b = reshape_in_layer_form(lv)

        epoch_cost = 0.0

        for batch in 1:n_batches
            i_start = (batch - 1) * batch_size + 1
            i_end   = min(batch * batch_size, n_samples)
            idx     = @view order[i_start:i_end]   # vue sur order, pas de copie

            Xb = @view Xtrain[idx, :]
            Yb = @view Ytrain[idx, :]

            # Forward pass — W, b passés directement
            Network.forward_pass!(buf, net, Xb, W, b)

            # Coût + dLF in-place
            J = _compute_cost_and_dLF!(dLF, adam.cost_function, Yb, buf.a_values[end])
            J += w2 * 0.5 * dot(θ, θ)

            # Backprop — W passé directement
            Network.backpropagation!(buf, net, Yb, dLF, W)
            @. buf.grad += w2 * θ

            # Mise à jour Adam fusionnée
            t  += 1
            @. m  = β1 * m + (1.0 - β1) * buf.grad
            @. v  = β2 * v + (1.0 - β2) * buf.grad^2
            α_t   = α * sqrt(1.0 - β2^t) / (1.0 - β1^t)
            @. θ -= α_t * m / (sqrt(v) + ε)

            epoch_cost += J
        end

        # W et b sont périmés après update θ — recalculés au début de la prochaine époque

        fplot[epoch] = epoch_cost / n_batches

        if epoch % 100 == 0 || epoch == 1
            @printf("Époque %4d / %d  |  coût = %.6g  |  %.1f s\n",
                    epoch, adam.max_epochs, fplot[epoch],
                    round(time() - start_time, digits=1))
        end
    end

    println("Entraînement terminé. Durée totale : $(round(time()-start_time, digits=2)) s")

    adam_final = AdamStruct(adam.cost_function, α, β1, β2, ε,
                             batch_size, adam.max_epochs, m, v, t, fplot)
    return adam_final, θ
end

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

function plot_cost_func(adam::AdamStruct)
    valid = adam.fplot[adam.fplot .> 0.0]
    isempty(valid) && (println("Aucune donnée de coût."); return)
    plt = plot(1:length(valid), valid;
               xlabel="Époques", ylabel="Coût", title="Convergence Adam",
               linewidth=1.8, legend=false, grid=true, yscale=:log10)
    display(plt)
    return plt
end

end