module Plot

using Plots
gr()

export plot_metrics,
       confusion_matrix_plot,
       roc_curve_plot,
       precision_recall_plot,
       f1_score

# Dossier de sortie pour tous les plots
const PLOT_DIR = "plots"

function ensure_plot_dir()
    isdir(PLOT_DIR) || mkdir(PLOT_DIR)
end

# =========================================================
# LOSS + ACCURACY
# =========================================================
function plot_metrics(train_losses, val_losses, train_accs, val_accs)
    ensure_plot_dir()

    p1 = plot(train_losses,
        label="Train Loss",
        xlabel="Epoch", ylabel="Loss",
        title="Loss Curve",
        linewidth=2)
    plot!(p1, val_losses, label="Validation Loss", linewidth=2)
    savefig(p1, joinpath(PLOT_DIR, "loss_curve.png"))
    println("  → Saved: $(PLOT_DIR)/loss_curve.png")

    p2 = plot(train_accs,
        label="Train Accuracy",
        xlabel="Epoch", ylabel="Accuracy",
        title="Accuracy Curve",
        linewidth=2)
    plot!(p2, val_accs, label="Validation Accuracy", linewidth=2)
    savefig(p2, joinpath(PLOT_DIR, "accuracy_curve.png"))
    println("  → Saved: $(PLOT_DIR)/accuracy_curve.png")
end

# =========================================================
# CONFUSION MATRIX
# =========================================================
function confusion_matrix_plot(y_true, y_pred)
    ensure_plot_dir()

    tp = sum((y_pred .== 1) .& (y_true .== 1))
    tn = sum((y_pred .== 0) .& (y_true .== 0))
    fp = sum((y_pred .== 1) .& (y_true .== 0))
    fn = sum((y_pred .== 0) .& (y_true .== 1))

    mat = [tn fp; fn tp]

    # annotations manuelles car annot=true n'est pas stable sur tous les backends
    ann = [(j, i, text(string(mat[i,j]), 14, :white, :center))
           for i in 1:2, j in 1:2]

    p = heatmap(mat,
        xticks=(1:2, ["Pred Neg", "Pred Pos"]),
        yticks=(1:2, ["True Neg", "True Pos"]),
        title="Confusion Matrix",
        c=:blues,
        annotations=ann,
        yflip=false)

    savefig(p, joinpath(PLOT_DIR, "confusion_matrix.png"))
    println("  → Saved: $(PLOT_DIR)/confusion_matrix.png")
    println("     TN=$tn  FP=$fp")
    println("     FN=$fn  TP=$tp")
end

# =========================================================
# F1 SCORE
# =========================================================
function f1_score(y_true, y_pred)

    tp = sum((y_pred .== 1) .& (y_true .== 1))
    fp = sum((y_pred .== 1) .& (y_true .== 0))
    fn = sum((y_pred .== 0) .& (y_true .== 1))

    precision = tp / (tp + fp + 1e-8)
    recall    = tp / (tp + fn + 1e-8)
    f1        = 2 * precision * recall / (precision + recall + 1e-8)

    println("Precision : ", round(precision, digits=4))
    println("Recall    : ", round(recall,    digits=4))
    println("F1-score  : ", round(f1,        digits=4))

    return f1
end

# =========================================================
# ROC CURVE
# =========================================================
function roc_curve_plot(y_true, y_scores)
    ensure_plot_dir()

    thresholds = collect(0:0.01:1)
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

    # AUC par méthode des trapèzes
    auc = -sum(diff(fpr_list) .* (tpr_list[1:end-1] .+ tpr_list[2:end]) ./ 2)

    p = plot(fpr_list, tpr_list,
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        title="ROC Curve (AUC = $(round(auc, digits=3)))",
        label="Model",
        legend=:bottomright,
        linewidth=2)
    plot!(p, [0,1], [0,1], linestyle=:dash, label="Random", color=:gray)

    savefig(p, joinpath(PLOT_DIR, "roc_curve.png"))
    println("  → Saved: $(PLOT_DIR)/roc_curve.png  (AUC = $(round(auc, digits=3)))")
end

# =========================================================
# PRECISION-RECALL CURVE
# =========================================================
function precision_recall_plot(y_true, y_scores)
    ensure_plot_dir()

    thresholds = collect(0:0.01:1)
    precisions = Float64[]
    recalls    = Float64[]

    for t in thresholds
        y_pred = y_scores .>= t
        tp = sum((y_pred .== 1) .& (y_true .== 1))
        fp = sum((y_pred .== 1) .& (y_true .== 0))
        fn = sum((y_pred .== 0) .& (y_true .== 1))
        push!(precisions, tp / (tp + fp + 1e-8))
        push!(recalls,    tp / (tp + fn + 1e-8))
    end

    p = plot(recalls, precisions,
        xlabel="Recall",
        ylabel="Precision",
        title="Precision-Recall Curve",
        label="Model",
        legend=:bottomleft,
        linewidth=2)

    savefig(p, joinpath(PLOT_DIR, "precision_recall.png"))
    println("  → Saved: $(PLOT_DIR)/precision_recall.png")
end

end