using Flux
using Statistics

include("Preprocess.jl")
include("Model.jl")
include("Train.jl")
include("Plot.jl")

using .Preprocess
using .Model
using .Train
using .Plot

# =========================================================
# PARAMÈTRES — tout ce qu'on veut changer est ici
# =========================================================

root_dir = joinpath(@__DIR__, "archive")

# ---------------------------------------------------------
# Choix du mode principal
#
# :custom → CNN from scratch (GroupNorm + GAP ou MLP)
# :resnet → ResNet18 pré-entraîné (Metalhead, RGB, 128×128)
# ---------------------------------------------------------
model_type = :custom # :custom ou :resnet

# ---------------------------------------------------------
# Paramètres mode :custom uniquement
# (ignorés si model_type = :resnet)
# ---------------------------------------------------------
img_size         = (128, 128)   # modifiable : (64,64) ou (128,128)
classifier       = :gap         # :gap ou :mlp
use_augmentation = false        # true si max_images < 2000

# ---------------------------------------------------------
# Paramètres communs aux deux modes
# ---------------------------------------------------------
batch_size = 32
epochs     = 30
lr         = 1e-3

# Smoke test   →  50
# Sanity check →  300
# Run complet  →  nothing
max_images = 3000

# =========================================================
# DATA — chargement selon le mode
# =========================================================

println("=== Chargement des données (mode $model_type) ===")

if model_type == :custom
    images, labels = load_dataset(root_dir;
                                  img_size=img_size,
                                  max_per_class=max_images)
else  # :resnet
    # img_size fixe 128×128 RGB — pas de paramètre à passer
    images, labels = load_dataset_resnet(root_dir;
                                         max_per_class=max_images)
end

train_data, val_data, test_data = split_dataset(images, labels)

# =========================================================
# MODEL
# =========================================================

println("\n=== Construction du modèle ===")

if model_type == :custom
    model = build_model(; model_type=:custom,
                          classifier=classifier,
                          img_size=img_size)
else
    model = build_model(; model_type=:resnet)
end

println(model)

# =========================================================
# TRAINING
# =========================================================

println("\n=== Entraînement ===")

# Pour ResNet en phase 1 (feature extraction) :
#   - LR plus petit car on affine seulement la tête
#   - Pas d'augmentation nécessaire (le backbone est robuste)
effective_lr  = model_type == :resnet ? lr * 0.1f0 : lr
effective_aug = model_type == :resnet ? false : use_augmentation

train_losses, val_losses, train_accs, val_accs =
    train_model(train_data, val_data, model;
                epochs=epochs,
                batch_size=batch_size,
                lr=effective_lr,
                patience=7,
                use_augmentation=effective_aug,
                model_type=model_type)

# =========================================================
# PLOTS TRAINING
# =========================================================

println("\n=== Visualisation ===")
plot_metrics(train_losses, val_losses, train_accs, val_accs)

# =========================================================
# ÉVALUATION FINALE
# =========================================================

println("\n=== Évaluation finale (test set) ===")

Flux.testmode!(model)

test_images, test_labels = test_data
y_scores = Float32[]

for x in test_images
    if model_type == :custom
        # grayscale : (H, W) → (W, H, 1, 1)
        x_in = reshape(Float32.(x), size(x,2), size(x,1), 1, 1)     #just to fit the correct type of data
    else
        # RGB : (H, W, 3) → (W, H, 3, 1)
        x_in = reshape(Float32.(x), size(x,2), size(x,1), 3, 1)
    end
    push!(y_scores, model(x_in)[1])
end

println("Scores sample : ", round.(y_scores[1:min(10,end)], digits=3))
println("Mean: $(round(mean(y_scores), digits=3)) | Min: $(round(minimum(y_scores), digits=3)) | Max: $(round(maximum(y_scores), digits=3))")

# Diagnostic par classe
pos_scores = y_scores[test_labels .== 1]
neg_scores = y_scores[test_labels .== 0]
println("\nPositifs (fissures)  → Mean: $(round(mean(pos_scores), digits=3))")
println("Négatifs (sans fiss) → Mean: $(round(mean(neg_scores), digits=3))")

sep = mean(pos_scores) - mean(neg_scores)
println("Séparation : $(round(sep, digits=4))")

if sep < -0.05
    println("⚠ Scores inversés → correction 1 - score appliquée")      #if mean pos < mean neg
    y_scores = 1f0 .- y_scores
end

# AUC
function compute_auc(y_true, y_scores)
    thresholds = collect(0.0:0.005:1.0)
    tpr_list = Float64[]
    fpr_list = Float64[]
    for t in thresholds
        y_pred = y_scores .>= t
        tp = sum((y_pred .== 1) .& (y_true .== 1))
        fp = sum((y_pred .== 1) .& (y_true .== 0))
        tn = sum((y_pred .== 0) .& (y_true .== 0))
        fn = sum((y_pred .== 0) .& (y_true .== 1))
        push!(tpr_list, tp / (tp + fn + 1e-8))        #recall or TPR
        push!(fpr_list, fp / (fp + tn + 1e-8))        #FPR
    end
    return -sum(diff(fpr_list) .* (tpr_list[1:end-1] .+ tpr_list[2:end]) ./ 2)
end

auc = compute_auc(test_labels, y_scores)
println("AUC : $(round(auc, digits=3))")

# Seuil optimal
function find_best_threshold(y_scores, y_true; thresholds=0.01:0.005:0.99)
    best_f1        = 0.0
    best_threshold = 0.5
    best_precision = 0.0
    best_recall    = 0.0
    for t in thresholds
        y_pred = y_scores .>= t
        tp = sum((y_pred .== 1) .& (y_true .== 1))
        fp = sum((y_pred .== 1) .& (y_true .== 0))
        fn = sum((y_pred .== 0) .& (y_true .== 1))
        prec = tp / (tp + fp + 1e-8)
        rec  = tp / (tp + fn + 1e-8)
        f1   = 2 * prec * rec / (prec + rec + 1e-8)
        if f1 > best_f1                 #we look for the best F1
            best_f1        = f1
            best_threshold = t
            best_precision = prec
            best_recall    = rec
        end
    end
    return best_threshold, best_precision, best_recall, best_f1
end

println("\n=== Recherche du seuil optimal ===")
best_threshold, best_precision, best_recall, best_f1 =
    find_best_threshold(y_scores, test_labels)

println("Seuil optimal : $(round(best_threshold, digits=3))")
println("Precision     : $(round(best_precision, digits=4))")
println("Recall        : $(round(best_recall,    digits=4))")
println("F1-score      : $(round(best_f1,        digits=4))")

y_pred_best = y_scores .>= best_threshold

f1_score(test_labels, y_pred_best)
confusion_matrix_plot(test_labels, y_pred_best)
roc_curve_plot(test_labels, y_scores)
precision_recall_plot(test_labels, y_scores)

println("\nDone.")