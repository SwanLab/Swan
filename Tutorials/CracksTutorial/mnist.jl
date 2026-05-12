# =========================================================
# mnist.jl
#
# Classification de chiffres manuscrits (0 à 9) avec un MLP
# Dataset : MNIST — 60 000 images 28×28 grayscale
#
# Objectif : valider que le pipeline Flux/Julia fonctionne
# sur un problème connu avant de conclure sur le projet fissures
#
# =========================================================

using Flux
using Flux: onehotbatch, onecold
using MLDatasets
using Statistics
using Random
using Optimisers
using Plots
gr()

Random.seed!(42)

const PLOT_DIR = joinpath(@__DIR__, "plots")

# =========================================================
# PARAMÈTRES
# =========================================================

batch_size = 64
epochs     = 20
lr         = 1e-3
patience   = 5

# =========================================================
# CHARGEMENT DES DONNÉES
#
# MNIST retourne :
#   x : Array Float32 de shape (28, 28, N)  — pixels dans [0,1]
#   y : Vector Int    de shape (N,)          — labels 0 à 9
# =========================================================

println("=== Chargement MNIST ===")

train_x, train_y = MNIST(:train)[:]
test_x,  test_y  = MNIST(:test)[:]

println("Train : $(size(train_x, 3)) images")
println("Test  : $(size(test_x,  3)) images")
println("Shape image : $(size(train_x)[1:2])")
println("Classes     : $(sort(unique(train_y)))")

# =========================================================
# PRÉPARATION DES DONNÉES
#
# Le MLP attend un vecteur 1D par image.
# On flatten chaque image 28×28 → vecteur de 784 valeurs.
#
# reshape(train_x, 784, :) :
#   (28, 28, 60000) → (784, 60000)
#   chaque colonne = une image aplatie
#
# onehotbatch : convertit les labels entiers en vecteurs binaires
#   label "3" → [0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
#   nécessaire pour la crossentropy catégorielle
# =========================================================

x_train = reshape(Float32.(train_x), 784, :)          # (784, 60000)
y_train = onehotbatch(train_y, 0:9)                    # (10, 60000)

x_test  = reshape(Float32.(test_x),  784, :)           # (784, 10000)
y_test  = onehotbatch(test_y,  0:9)                    # (10, 10000)

println("\nShape x_train : $(size(x_train))")
println("Shape y_train : $(size(y_train))")

# =========================================================
# MODÈLE — MLP simple
#
# 784 → 256 → 128 → 10
#
# Différences vs projet fissures :
#   - Pas de couches Conv (images petites et bien centrées)
#   - softmax en sortie au lieu de sigmoid (10 classes)
#   - crossentropy catégorielle au lieu de binaire
# =========================================================

println("\n=== Construction du modèle ===")

model = Chain(
    Dense(784, 256, relu),
    Dropout(0.3),
    Dense(256, 128, relu),
    Dropout(0.2),
    Dense(128, 10),
    softmax                   # probabilité pour chacun des 10 chiffres
)

println(model)
nb_params = sum(length, Flux.params(model))
println("Nombre de paramètres : $nb_params")

# =========================================================
# FONCTIONS UTILITAIRES
# =========================================================

# Accuracy : compare la classe prédite (argmax) au label réel
function accuracy(model, x, y)
    ŷ     = model(x)
    preds = onecold(ŷ, 0:9)    # argmax → chiffre prédit
    trues = onecold(y, 0:9)    # argmax → chiffre réel
    return mean(preds .== trues)
end

# Loss : crossentropy catégorielle
# -Σ y_i × log(ŷ_i) sur les 10 classes
# Seulement la classe correcte contribue (les autres ont y_i=0)
loss_fn(m, x, y) = Flux.crossentropy(m(x), y)

# =========================================================
# ENTRAÎNEMENT
# =========================================================

println("\n=== Entraînement ===")

opt = Optimisers.Adam(lr)
st  = Optimisers.setup(opt, model)

N = size(x_train, 2)   # 60 000

train_losses = Float32[]
val_losses   = Float32[]
train_accs   = Float32[]
val_accs     = Float32[]

best_val_loss    = Inf32
patience_counter = 0

for epoch in 1:epochs

    Flux.trainmode!(model)

    idxs       = shuffle(1:N)
    epoch_loss = 0f0
    nb_batches = 0

    for i in 1:batch_size:N
        batch_idxs = idxs[i:min(i+batch_size-1, N)]

        x_batch = x_train[:, batch_idxs]   # (784, batch_size)
        y_batch = y_train[:, batch_idxs]   # (10,  batch_size)

        gs = gradient(m -> loss_fn(m, x_batch, y_batch), model)[1]
        st, model = Optimisers.update(st, model, gs)

        epoch_loss += Flux.crossentropy(model(x_batch), y_batch)
        nb_batches += 1
    end

    train_loss = epoch_loss / nb_batches

    # Validation sur le test set complet
    Flux.testmode!(model)
    val_loss = Flux.crossentropy(model(x_test), y_test)
    tr_acc   = accuracy(model, x_train, y_train)
    va_acc   = accuracy(model, x_test,  y_test)

    push!(train_losses, train_loss)
    push!(val_losses,   val_loss)
    push!(train_accs,   tr_acc)
    push!(val_accs,     va_acc)

    println("Epoch $epoch | Train Loss: $(round(train_loss, digits=4)) Acc: $(round(tr_acc*100, digits=1))% | Val Loss: $(round(val_loss, digits=4)) Acc: $(round(va_acc*100, digits=1))%")

    if val_loss < best_val_loss
        best_val_loss    = val_loss
        patience_counter = 0
    else
        patience_counter += 1
    end

    if patience_counter >= patience
        println("Early stopping à l'epoch $epoch")
        break
    end
end

# =========================================================
# PLOTS
# =========================================================

isdir(PLOT_DIR) || mkdir(PLOT_DIR)

p1 = plot(train_losses, label="Train Loss", linewidth=2,
          xlabel="Epoch", ylabel="Loss", title="Loss Curve")
plot!(p1, val_losses, label="Val Loss", linewidth=2)
savefig(p1, joinpath(PLOT_DIR, "mnist_loss.png"))
println("\n→ Saved: mnist_loss.png")

p2 = plot(train_accs .* 100, label="Train Acc", linewidth=2,
          xlabel="Epoch", ylabel="Accuracy (%)", title="Accuracy Curve")
plot!(p2, val_accs .* 100, label="Val Acc", linewidth=2)
savefig(p2, joinpath(PLOT_DIR, "mnist_accuracy.png"))
println("→ Saved: mnist_accuracy.png")

# =========================================================
# ÉVALUATION FINALE
# =========================================================

println("\n=== Évaluation finale ===")

Flux.testmode!(model)

final_acc = accuracy(model, x_test, y_test)
println("Accuracy test set : $(round(final_acc*100, digits=2))%")

# Matrice de confusion — montre quels chiffres sont confondus
println("\nMatrice de confusion (lignes=réel, colonnes=prédit) :")

ŷ_classes = onecold(model(x_test), 0:9)
y_classes  = onecold(y_test, 0:9)

conf_matrix = zeros(Int, 10, 10)
for (pred, true_) in zip(ŷ_classes, y_classes)
    conf_matrix[true_+1, pred+1] += 1
end

# Affichage lisible
print("     ")
for j in 0:9
    print("  $j  ")
end
println()
for i in 1:10
    print("  $(i-1)  ")
    for j in 1:10
        val = conf_matrix[i,j]
        if i == j
            print("[$val]")   # diagonal = bonnes prédictions
        else
            print(" $val ")
        end
    end
    println()
end

p3 = heatmap(conf_matrix,
    xticks=(1:10, string.(0:9)),
    yticks=(1:10, string.(0:9)),
    xlabel="Prédit", ylabel="Réel",
    title="Matrice de confusion MNIST",
    c=:blues)
savefig(p3, joinpath(PLOT_DIR, "mnist_confusion.png"))
println("\n→ Saved: mnist_confusion.png")

println("\nDone.")