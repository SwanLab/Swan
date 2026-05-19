# =========================================================
# mnist_manual.jl
#
# Classification MNIST avec MLP entièrement fait à la main.
# Forward pass, backward pass et Adam implémentés manuellement,
# calqués sur l'architecture LearnableVariables / Network / Optimizer
# du projet de régression.
#
# Seul Flux/MLDatasets est utilisé pour le chargement des images.
#
# Architecture : 784 → 256 → 128 → 10
# Activation   : ReLU (couches cachées), Softmax (sortie)
# Loss         : Cross-entropie catégorielle
# Optimiseur   : Adam
# =========================================================

using MLDatasets
using LinearAlgebra
using Random
using Distributions
using Statistics
using Printf
using Plots
gr()

Random.seed!(42)

const PLOT_DIR = joinpath(@__DIR__, "plots")
isdir(PLOT_DIR) || mkdir(PLOT_DIR)

# =========================================================
# HYPERPARAMÈTRES
# =========================================================

const HIDDEN_LAYERS = [256, 128]
const BATCH_SIZE    = 128
const MAX_EPOCHS    = 50
const LR            = 1e-3
const β1            = 0.9
const β2            = 0.999
const ε_adam        = 1e-8
const λ             = 1e-4      # régularisation L2
const PATIENCE      = 8         # early stopping

# =========================================================
# CHARGEMENT ET PRÉPARATION DES DONNÉES
#
# MNIST : images 28×28, labels 0-9
# On flatten chaque image → vecteur 784
# Les données sont en Float64 pour cohérence avec le projet regression
# =========================================================

println("=== Chargement MNIST ===")

train_x_raw, train_y_raw = MNIST(:train)[:]
test_x_raw,  test_y_raw  = MNIST(:test)[:]

# Flatten + transposer : on veut (N, 784) — lignes = exemples, comme dans le projet régression
# Le projet régression utilise (N, n_features) pour Xtrain
Xtrain = Float64.(reshape(train_x_raw, 784, :)')   # (60000, 784)
Xtest  = Float64.(reshape(test_x_raw,  784, :)')   # (10000, 784)

# One-hot : (N, 10) — même convention que le projet régression (N, n_labels)
function one_hot(labels::Vector, n_classes::Int)
    N = length(labels)
    Y = zeros(Float64, N, n_classes)
    for (i, l) in enumerate(labels)
        Y[i, l + 1] = 1.0   # labels MNIST démarrent à 0
    end
    return Y
end

Ytrain = one_hot(train_y_raw, 10)   # (60000, 10)
Ytest  = one_hot(test_y_raw,  10)   # (10000, 10)

println("Train : $(size(Xtrain, 1)) images — shape $(size(Xtrain))")
println("Test  : $(size(Xtest,  1)) images — shape $(size(Xtest))")
println("Labels train : shape $(size(Ytrain))")

# =========================================================
# LEARNABLE VARIABLES
#
# Reprise exacte de la logique du module LearnableVariables :
# - Tous les paramètres (W, b) aplatis dans un vecteur θ
# - Initialisation Xavier uniform pour W, 0.1 pour b (hidden), 1/n pour sortie
# - reshape_in_layer_form retourne des vues (@view) sans copie
# =========================================================

struct LearnableVars
    neurons_per_layer :: Vector{Int}
    n_layers          :: Int
    thetavec          :: Vector{Float64}
end

function init_learnable_variables(neurons_per_layer::Vector{Int})
    npl      = neurons_per_layer
    n_layers = length(npl)
    total    = sum(npl[i-1] * npl[i] + npl[i] for i in 2:n_layers)
    θ        = Vector{Float64}(undef, total)

    offset = 1
    for i in 2:n_layers
        in_dim, out_dim = npl[i-1], npl[i]
        nW = in_dim * out_dim

        # Xavier uniform
        u = sqrt(6.0 / (in_dim + out_dim))
        θ[offset : offset + nW - 1] .= rand(Uniform(-u, u), nW)

        # Biais
        b_val = (i != n_layers) ? 0.1 : 1.0 / out_dim
        θ[offset + nW : offset + nW + out_dim - 1] .= b_val

        offset += nW + out_dim
    end

    return LearnableVars(npl, n_layers, θ)
end

"""
    reshape_in_layer_form(lv)

Retourne (W, b) comme vecteurs de matrices/vecteurs.
W[i] : matrice (npl[i], npl[i+1]) — convention (in, out) pour mul!(Y, X, W)
b[i] : vecteur (npl[i+1],)
Utilise @view — zéro copie.
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
        nW       = in_size * out_size

        W[i] = reshape(@view(θ[offset : offset + nW - 1]), in_size, out_size)
        b[i] = @view(θ[offset + nW : offset + nW + out_size - 1])

        offset += nW + out_size
    end
    return W, b
end

function update_thetavec(lv::LearnableVars, θ_new::Vector{Float64})
    return LearnableVars(lv.neurons_per_layer, lv.n_layers, θ_new)
end

# =========================================================
# ACTIVATIONS
#
# ReLU pour les couches cachées — identique au projet régression.
# Softmax pour la sortie — remplace linear/sigmoid.
#   softmax(z)_i = exp(z_i) / Σ exp(z_j)
#   Stabilité numérique : on soustrait max(z) avant exp (log-sum-exp trick)
# =========================================================

function relu!(a::AbstractMatrix{Float64})
    @. a = max(a, 0.0)
end

function relu_derivative!(out::AbstractMatrix{Float64}, a::AbstractMatrix{Float64})
    @. out = Float64(a > 0.0)
end

"""
    softmax!(a)

Softmax in-place, ligne par ligne (chaque ligne = un exemple).
Soustrait le max par ligne pour la stabilité numérique.
"""
function softmax!(a::AbstractMatrix{Float64})
    @inbounds for i in axes(a, 1)
        row_max = -Inf
        for j in axes(a, 2)
            row_max = max(row_max, a[i, j])
        end
        s = 0.0
        for j in axes(a, 2)
            a[i, j] = exp(a[i, j] - row_max)
            s += a[i, j]
        end
        for j in axes(a, 2)
            a[i, j] /= s
        end
    end
end

# =========================================================
# NETWORK BUFFERS
#
# Reprise exacte de NetworkBuffers du projet régression.
# Tous les buffers sont alloués une seule fois et réutilisés.
# =========================================================

struct NetworkBuffers
    a_values :: Vector{Matrix{Float64}}   # activations par couche
    deltag   :: Vector{Matrix{Float64}}   # δ (erreur rétropropagée)
    dcW      :: Vector{Matrix{Float64}}   # gradient ∂C/∂W par couche
    dcB      :: Vector{Vector{Float64}}   # gradient ∂C/∂b par couche
    grad     :: Vector{Float64}           # gradient aplati (même taille que θ)
    tmp      :: Vector{Matrix{Float64}}   # buffer intermédiaire backprop
end

function init_buffers(npl::Vector{Int}, n_layers::Int, batch_size::Int)
    a_values = Vector{Matrix{Float64}}(undef, n_layers)
    a_values[1] = zeros(batch_size, npl[1])
    for i in 2:n_layers
        a_values[i] = zeros(batch_size, npl[i])
    end

    deltag = Vector{Matrix{Float64}}(undef, n_layers)
    for k in 2:n_layers
        deltag[k] = zeros(batch_size, npl[k])
    end

    dcW  = [zeros(npl[i], npl[i+1]) for i in 1:n_layers-1]
    dcB  = [zeros(npl[i+1])          for i in 1:n_layers-1]
    grad = zeros(sum(npl[i]*npl[i+1] + npl[i+1] for i in 1:n_layers-1))

    tmp = Vector{Matrix{Float64}}(undef, n_layers)
    for k in 2:n_layers
        tmp[k] = zeros(batch_size, npl[k])
    end

    return NetworkBuffers(a_values, deltag, dcW, dcB, grad, tmp)
end

# =========================================================
# FORWARD PASS
#
# Identique au projet régression, sauf activation de sortie :
#   - couches cachées : ReLU (in-place)
#   - couche de sortie : Softmax (in-place)
#
# Convention : X est (batch_size, n_features) — lignes = exemples
# mul!(Y, X, W) calcule X @ W sans allocation
# b[i]' broadcaste le vecteur biais sur toutes les lignes
# =========================================================

function forward_pass!(buf::NetworkBuffers, npl::Vector{Int}, n_layers::Int,
                        X::AbstractMatrix{Float64},
                        W::Vector{<:AbstractMatrix{Float64}},
                        b::Vector{<:AbstractVector{Float64}})
    buf.a_values[1] = X   # vue, pas de copie

    for i in 2:n_layers
        mul!(buf.a_values[i], buf.a_values[i-1], W[i-1])
        buf.a_values[i] .+= b[i-1]'

        if i == n_layers
            softmax!(buf.a_values[i])
        else
            relu!(buf.a_values[i])
        end
    end
end

# =========================================================
# LOSS : CROSS-ENTROPIE CATÉGORIELLE + GRADIENT
#
# J = -1/m Σ_i Σ_k Y[i,k] * log(Ŷ[i,k])
#
# Gradient de la loss par rapport à ŷ AVANT softmax
# (dérivée combinée softmax + cross-entropie) :
#   dL/dz = (ŷ - y) / m
#
# Ce résultat simplifié découle du fait que pour softmax + cross-entropie,
# le Jacobien de softmax et la dérivée de la log-loss se combinent en (ŷ-y).
# Dans le code, on passe ce gradient directement comme deltag[nL],
# ce qui court-circuite la multiplication par la dérivée de softmax
# (qui serait mal définie avec la convention scalaire du projet régression).
# =========================================================

function crossentropy_and_gradient!(dLF::Matrix{Float64},
                                     Y::AbstractMatrix{Float64},
                                     Ŷ::AbstractMatrix{Float64})
    m = size(Y, 1)
    J = 0.0
    @inbounds for j in axes(Y, 2), i in axes(Y, 1)
        yhat = max(Ŷ[i, j], 1e-15)   # stabilité numérique
        J   -= Y[i, j] * log(yhat)
    end
    J /= m

    # Gradient simplifié softmax + cross-entropie : (ŷ - y) / m
    @. dLF = (Ŷ - Y) / m

    return J
end

# =========================================================
# BACKWARD PASS (BACKPROPAGATION)
#
# Reprend exactement la structure de backpropagation! du projet régression.
#
# Différence clé pour la couche de sortie (softmax) :
#   Dans le projet régression, on calcule d'abord la dérivée de l'activation,
#   puis on la multiplie par dLF.
#   Ici, pour softmax + cross-entropie, le gradient combiné est directement
#   dLF = (ŷ - y)/m (déjà calculé dans crossentropy_and_gradient!).
#   On affecte donc directement deltag[nL] = dLF sans passer par
#   _activation_derivative, ce qui est mathématiquement équivalent.
#
# Pour les couches cachées (ReLU), la logique est identique au projet :
#   δ[k] = (δ[k+1] @ W[k]') ⊙ relu'(a[k])
# =========================================================

function backpropagation!(buf::NetworkBuffers, npl::Vector{Int}, n_layers::Int,
                           Y::AbstractMatrix{Float64},
                           dLF::Matrix{Float64},
                           W::Vector{<:AbstractMatrix{Float64}})
    m  = size(Y, 1)
    nL = n_layers

    for k in reverse(2:nL)
        if k == nL
            # Couche de sortie : gradient combiné softmax+crossentropy déjà dans dLF
            @. buf.deltag[k] = dLF
        else
            # Couches cachées : ReLU
            relu_derivative!(buf.deltag[k], buf.a_values[k])
            mul!(buf.tmp[k], buf.deltag[k+1], W[k]')
            buf.deltag[k] .*= buf.tmp[k]
        end

        # Gradient par rapport aux poids : X_k' @ δ_k+1 / m
        mul!(buf.dcW[k-1], buf.a_values[k-1]', buf.deltag[k])
        buf.dcW[k-1] ./= m

        # Gradient par rapport aux biais : somme sur le batch / m
        sum!(reshape(buf.dcB[k-1], 1, :), buf.deltag[k])
        buf.dcB[k-1] ./= m
    end

    # Aplatissement in-place dans buf.grad (même structure que le projet régression)
    offset = 1
    for i in 1:nL-1
        nW = npl[i] * npl[i+1]
        buf.grad[offset : offset + nW - 1]                  .= vec(buf.dcW[i])
        buf.grad[offset + nW : offset + nW + npl[i+1] - 1] .= buf.dcB[i]
        offset += nW + npl[i+1]
    end
end

# =========================================================
# PRÉCISION (ACCURACY)
#
# Comparaison argmax(ŷ) vs argmax(y) sur tout un dataset.
# Appelé hors boucle d'entraînement — allocation autorisée.
# =========================================================

function accuracy(lv::LearnableVars, X::Matrix{Float64}, Y::Matrix{Float64})
    W, b = reshape_in_layer_form(lv)
    npl  = lv.neurons_per_layer
    nL   = lv.n_layers
    a    = X
    for i in 2:nL
        z = a * W[i-1] .+ b[i-1]'
        if i == nL
            # softmax
            z .-= maximum(z, dims=2)
            z   = exp.(z)
            z ./= sum(z, dims=2)
            a   = z
        else
            a = max.(z, 0.0)   # ReLU
        end
    end
    preds  = argmax.(eachrow(a)) .- 1    # index julia 1-based → label 0-based
    truths = argmax.(eachrow(Y)) .- 1
    return mean(preds .== truths)
end

# =========================================================
# ADAM — BOUCLE D'ENTRAÎNEMENT
#
# Reprend _training_loop du projet régression avec :
# - Early stopping sur la loss de validation
# - Calcul de l'accuracy train/val chaque époque
# - Régularisation L2 sur θ (grad += λ*θ)
#
# Convention :
#   α_t = α * sqrt(1 - β2^t) / (1 - β1^t)
#   m   = β1*m + (1-β1)*grad
#   v   = β2*v + (1-β2)*grad²
#   θ  -= α_t * m / (sqrt(v) + ε)
# =========================================================

function train_adam!(lv::LearnableVars,
                     Xtrain::Matrix{Float64}, Ytrain::Matrix{Float64},
                     Xtest::Matrix{Float64},  Ytest::Matrix{Float64};
                     batch_size=128, max_epochs=50, lr=1e-3,
                     β1=0.9, β2=0.999, ε=1e-8, λ=1e-4, patience=8)

    npl      = lv.neurons_per_layer
    n_layers = lv.n_layers
    θ        = copy(lv.thetavec)
    n_params = length(θ)
    n_train  = size(Xtrain, 1)
    n_inputs  = size(Xtrain, 2)
    n_outputs = size(Ytrain, 2)

    # États Adam
    m_adam = zeros(n_params)
    v_adam = zeros(n_params)
    t_adam = 0

    # Buffers pré-alloués (même logique que init_buffers du projet régression)
    buf   = init_buffers(npl, n_layers, batch_size)
    dLF   = zeros(batch_size, n_outputs)
    order = Vector{Int}(undef, n_train)
    Xb_buf = Matrix{Float64}(undef, batch_size, n_inputs)
    Yb_buf = Matrix{Float64}(undef, batch_size, n_outputs)

    # Historiques
    train_losses = Float64[]
    val_losses   = Float64[]
    train_accs   = Float64[]
    val_accs     = Float64[]

    best_val_loss    = Inf
    patience_counter = 0
    start_time       = time()

    for epoch in 1:max_epochs
        randperm!(order)
        n_batches  = max(1, fld(n_train, batch_size))
        epoch_loss = 0.0

        # Extraction des vues W, b depuis θ courant
        lv_cur   = update_thetavec(lv, θ)
        W, b     = reshape_in_layer_form(lv_cur)

        for batch in 1:n_batches
            i_start   = (batch - 1) * batch_size + 1
            i_end     = min(batch * batch_size, n_train)
            actual_bs = i_end - i_start + 1

            # Copie dans buffers contigus (comme dans _training_loop)
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

            # Forward
            forward_pass!(buf, npl, n_layers, Xb, W, b)

            # Loss + gradient de la loss (cross-entropie)
            dLF_view = @view dLF[1:actual_bs, :]
            J = crossentropy_and_gradient!(dLF_view, Yb, buf.a_values[end])
            J += 0.5 * λ * dot(θ, θ)

            # Backward
            backpropagation!(buf, npl, n_layers, Yb, dLF_view, W)

            # Régularisation L2 : grad += λ*θ
            @. buf.grad += λ * θ

            # Mise à jour Adam
            t_adam += 1
            α_t = lr * sqrt(1.0 - β2^t_adam) / (1.0 - β1^t_adam)
            @. m_adam = β1 * m_adam + (1.0 - β1) * buf.grad
            @. v_adam = β2 * v_adam + (1.0 - β2) * buf.grad^2
            @. θ -= α_t * m_adam / (sqrt(v_adam) + ε)

            # W et b sont des vues sur l'ancien θ — on les recalcule
            lv_cur = update_thetavec(lv, θ)
            W, b   = reshape_in_layer_form(lv_cur)

            epoch_loss += J
        end

        avg_loss = epoch_loss / n_batches

        # Validation (hors boucle de batch — allocation OK)
        lv_eval  = update_thetavec(lv, θ)
        W_e, b_e = reshape_in_layer_form(lv_eval)

        buf_val = init_buffers(npl, n_layers, size(Xtest, 1))
        dLF_val = zeros(size(Xtest, 1), n_outputs)
        forward_pass!(buf_val, npl, n_layers, Xtest, W_e, b_e)
        val_J = crossentropy_and_gradient!(dLF_val, Ytest, buf_val.a_values[end])

        tr_acc  = accuracy(lv_eval, Xtrain, Ytrain)
        va_acc  = accuracy(lv_eval, Xtest,  Ytest)

        push!(train_losses, avg_loss)
        push!(val_losses,   val_J)
        push!(train_accs,   tr_acc)
        push!(val_accs,     va_acc)

        @printf("Epoch %3d/%d | Train Loss: %.4f  Acc: %.1f%%  | Val Loss: %.4f  Acc: %.1f%%  | %.1f s\n",
                epoch, max_epochs, avg_loss, tr_acc*100, val_J, va_acc*100,
                time() - start_time)

        # Early stopping
        if val_J < best_val_loss
            best_val_loss    = val_J
            patience_counter = 0
        else
            patience_counter += 1
        end

        if patience_counter >= patience
            println("Early stopping à l'epoch $epoch (patience=$patience)")
            break
        end
    end

    println("\nEntraînement terminé. Durée : $(round(time()-start_time, digits=1)) s")
    return θ, train_losses, val_losses, train_accs, val_accs
end

# =========================================================
# INITIALISATION ET ENTRAÎNEMENT
# =========================================================

println("\n=== Initialisation du réseau ===")

neurons_per_layer = [784; HIDDEN_LAYERS; 10]
println("Architecture : $(join(neurons_per_layer, " → "))")

lv = init_learnable_variables(neurons_per_layer)
n_params = length(lv.thetavec)
println("Nombre de paramètres : $n_params")

println("\n=== Entraînement (Adam) ===\n")

θ_final, train_losses, val_losses, train_accs, val_accs = train_adam!(
    lv, Xtrain, Ytrain, Xtest, Ytest;
    batch_size = BATCH_SIZE,
    max_epochs = MAX_EPOCHS,
    lr         = LR,
    β1         = β1,
    β2         = β2,
    ε          = ε_adam,
    λ          = λ,
    patience   = PATIENCE
)

# =========================================================
# PLOTS
# =========================================================

p1 = plot(train_losses, label="Train Loss", linewidth=2,
          xlabel="Epoch", ylabel="Loss (cross-entropie)", title="Loss Curve")
plot!(p1, val_losses, label="Val Loss", linewidth=2)
savefig(p1, joinpath(PLOT_DIR, "mnist_manual_loss.png"))
println("\n→ Saved: mnist_manual_loss.png")

p2 = plot(train_accs .* 100, label="Train Acc", linewidth=2,
          xlabel="Epoch", ylabel="Accuracy (%)", title="Accuracy Curve")
plot!(p2, val_accs .* 100, label="Val Acc", linewidth=2)
savefig(p2, joinpath(PLOT_DIR, "mnist_manual_accuracy.png"))
println("→ Saved: mnist_manual_accuracy.png")

# =========================================================
# ÉVALUATION FINALE
# =========================================================

println("\n=== Évaluation finale ===")

lv_final = update_thetavec(lv, θ_final)
final_acc = accuracy(lv_final, Xtest, Ytest)
println("Accuracy test set : $(round(final_acc*100, digits=2))%")

# Matrice de confusion
println("\nMatrice de confusion (lignes=réel, colonnes=prédit) :")

W_f, b_f = reshape_in_layer_form(lv_final)
buf_conf  = init_buffers(lv_final.neurons_per_layer, lv_final.n_layers, size(Xtest, 1))
forward_pass!(buf_conf, lv_final.neurons_per_layer, lv_final.n_layers, Xtest, W_f, b_f)

ŷ_classes = [argmax(buf_conf.a_values[end][i, :]) - 1 for i in 1:size(Xtest, 1)]
y_classes  = [argmax(Ytest[i, :]) - 1               for i in 1:size(Xtest, 1)]

conf_matrix = zeros(Int, 10, 10)
for (pred, true_) in zip(ŷ_classes, y_classes)
    conf_matrix[true_+1, pred+1] += 1
end

print("     ")
for j in 0:9; print("  $j  "); end
println()
for i in 1:10
    print("  $(i-1)  ")
    for j in 1:10
        val = conf_matrix[i, j]
        print(i == j ? "[$val]" : " $val ")
    end
    println()
end

p3 = heatmap(conf_matrix,
    xticks=(1:10, string.(0:9)),
    yticks=(1:10, string.(0:9)),
    xlabel="Prédit", ylabel="Réel",
    title="Matrice de confusion MNIST (manuel)",
    c=:blues)
savefig(p3, joinpath(PLOT_DIR, "mnist_manual_confusion.png"))
println("\n→ Saved: mnist_manual_confusion.png")

println("\nDone.")