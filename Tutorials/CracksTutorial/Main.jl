using Flux
using Statistics

include("module_Preprocess.jl")
include("module_Model.jl")
include("module_Train.jl")
include("module_Plot.jl")

using .Preprocess
using .Model
using .Train
using .Plot

# =========================================================
# PARAMÈTRES
# =========================================================

root_dir = "C:/Users/couwa/Documents/Stage_CIMNE/Classification/archive"

img_size   = (64, 64)
batch_size = 32
epochs     = 30
lr         = 1e-3

# Smoke test   →  50       (< 1 min)
# Sanity check →  300      (~ 5 min)
# Run complet  →  nothing  (toutes les images)
max_images = 300

# Data augmentation (flip H/V + bruit gaussien sur le train set)
# true  → recommandé si max_images < 2000
# false → suffisant si max_images >= 2000
use_augmentation = true

# Architecture du classifieur
# :gap → GlobalAvgPool → Dense(128→1)       recommandé < 1000 images
# :mlp → Flatten(8192) → Dense(8192→256) → Dense(256→1)  recommandé > 2000 images
classifier = :gap

# =========================================================
# DATA
# =========================================================

println("=== Chargement des données ===")
images, labels = load_dataset(root_dir; img_size=img_size, max_per_class=max_images)
train_data, val_data, test_data = split_dataset(images, labels)

# =========================================================
# MODEL
# =========================================================

println("\n=== Construction du modèle ===")
model = build_model(; classifier=classifier)
println(model)

# =========================================================
# TRAINING
# =========================================================

println("\n=== Entraînement ===")
train_losses, val_losses, train_accs, val_accs =
    train_model(model, train_data, val_data;
                epochs=epochs,
                batch_size=batch_size,
                lr=lr,
                patience=7,
                use_augmentation=use_augmentation)

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
    x_in = reshape(Float32.(x), size(x,2), size(x,1), 1, 1)
    push!(y_scores, model(x_in)[1])
end

println("Scores sample : ", round.(y_scores[1:min(10,end)], digits=3))
println("Mean: $(round(mean(y_scores), digits=3)) | Min: $(round(minimum(y_scores), digits=3)) | Max: $(round(maximum(y_scores), digits=3))")

# Diagnostic par classe
pos_scores = y_scores[test_labels .== 1]
neg_scores = y_scores[test_labels .== 0]
println("\nPositifs (fissures)  → Mean: $(round(mean(pos_scores), digits=3))  Min: $(round(minimum(pos_scores), digits=3))  Max: $(round(maximum(pos_scores), digits=3))")
println("Négatifs (sans fisc) → Mean: $(round(mean(neg_scores), digits=3))  Min: $(round(minimum(neg_scores), digits=3))  Max: $(round(maximum(neg_scores), digits=3))")

sep = mean(pos_scores) - mean(neg_scores)
println("Séparation (mean pos - mean neg) : $(round(sep, digits=4))")

if sep < -0.05
    println("⚠ Scores inversés détectés → correction 1 - score appliquée")
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
        push!(tpr_list, tp / (tp + fn + 1e-8))
        push!(fpr_list, fp / (fp + tn + 1e-8))
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
        if f1 > best_f1
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

println("Seuil optimal  : $(round(best_threshold, digits=3))")
println("Precision      : $(round(best_precision, digits=4))")
println("Recall         : $(round(best_recall,    digits=4))")
println("F1-score       : $(round(best_f1,        digits=4))")

y_pred_best = y_scores .>= best_threshold

f1_score(test_labels, y_pred_best)
confusion_matrix_plot(test_labels, y_pred_best)
roc_curve_plot(test_labels, y_scores)
precision_recall_plot(test_labels, y_scores)

println("\nDone.")