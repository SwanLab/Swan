module LossFunctional

export LossFunctional, init_lossfunctional, compute_loss_and_gradient, compute_stochastic_loss_and_gradient, get_test_error

using LinearAlgebra
using Network
using LearnableVariables
using Data
using Random

"""
    LossFunctional

A container for data, cost type, network, and learnable variables.
Immutable.
"""
struct LossFunctionalStruct
    cost_type::String
    data::DataStruct
    network::Net
end

"""
    init_lossfunctional(params)

Constructs a LossFunctional from parameter dictionary.
"""
function init_lossfunctional(params::Dict{String,Any})
    return LossFunctionalStruct(
        params["costType"],
        params["data"],
        params["network"]
    )
end

"""
    compute_loss_and_gradient(lf, θ)

Full-batch cost and gradient.
This includes auxiliar function computeGradient()
"""
function compute_loss_and_gradient(lf::LossFunctionalStruct, θ::Vector{Float64})
    lf.network.learnable_variables = set_theta(lf.network.learnable_variables, θ)
    X = lf.data.Xtrain
    Y = lf.data.Ytrain

    a_vals = forward_pass(lf.network, X)
    y_out = a_vals[end]
    j, dLF = compute_cost_and_dLF(lf, Y, y_out)
    grad = backpropagation(lf.network, Y, dLF, a_vals)

    return j, grad
end

"""
    compute_stochastic_loss_and_gradient(lf, x, i_batch, order, move_batch)

Mini-batch cost and gradient with batch update logic.
Returns (loss, grad, is_depleted, next_batch_id).
"""
function compute_stochastic_loss_and_gradient(lf::LossFunctionalStruct, θ::Vector{Float64}, i_batch::Int, order::Vector{Int}, move_batch::Bool)
    lf.network.learnable_variables = set_theta(lf.network.learnable_variables, θ)

    X = lf.data.Xtrain
    Y = lf.data.Ytrain
    batch_size = compute_batch_size(X)
    n_batches = fld(size(X, 1), batch_size)

    i_start = (i_batch - 1) * batch_size + 1
    i_end = min(i_batch * batch_size, size(X, 1))
    idx = order[i_start:i_end]

    Xb = X[idx, :]
    Yb = Y[idx, :]

    z_vals, a_vals = forward_pass(lf.network, Xb)
    y_out = a_vals[end]
    j, dLF = compute_cost_and_dLF(lf, Yb, y_out)
    grad = backpropagation(lf.network, Yb, dLF, a_vals)

    is_depleted = move_batch && i_batch == n_batches
    #next_batch = move_batch ? (i_batch % n_batches) + 1 : i_batch

    return j, grad, is_depleted #, next_batch
end

"""
    get_test_error(lf)

Returns misclassification rate on test set.
"""
function get_test_error(lf::LossFunctionalStruct)
    H = compute_last_H(lf.network, lf.data.Xtest)
    Ypred = map(row -> findmax(row)[2], eachrow(H))
    Ytarget = map(row -> findmax(row)[2], eachrow(lf.data.Ytest))
    return sum(Ypred .!= Ytarget) / length(Ytarget)
end

"""
    compute_cost_and_dLF(lf, Y, Ŷ)

Returns (loss, derivative w.r.t. final layer).
"""
function compute_cost_and_dLF(lf::LossFunctional, Y::Matrix{Float64}, Ŷ::Matrix{Float64})
    yp = Ŷ .- 1e-11  # for stability
    if lf.cost_type == "-loglikelihood"
        c = sum((1 .- Y) .* (-log.(1 .- yp)) .+ Y .* (-log.(yp)), dims=2)
        J = mean(c)
        dLF = (yp .- Y) ./ (yp .* (1 .- yp))
    elseif lf.cost_type == "L2"
        c = (Ŷ .- Y).^2
        J = sqrt(sum(c))
        dLF = Ŷ .- Y
    else
        error("Invalid loss type: $(lf.cost_type)")
    end
    return J, dLF
end

"""
    compute_batch_size(X)

Batch size policy (cap at 200).
"""
compute_batch_size(X::Matrix) = min(size(X, 1), 200)

end
