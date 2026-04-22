module Plot

using Plots

export plot_metrics,
       confusion_matrix_plot,
       roc_curve_plot,
       precision_recall_plot,
       f1_score

# =========================================================
# 🔹 LOSS + ACCURACY
# =========================================================
function plot_metrics(train_losses, val_losses,
                      train_accs, val_accs)

    p1 = plot(train_losses,
        label="Train Loss",
        xlabel="Epoch",
        ylabel="Loss",
        title="Loss Curve")

    plot!(p1, val_losses, label="Validation Loss")

    p2 = plot(train_accs,
        label="Train Accuracy",
        xlabel="Epoch",
        ylabel="Accuracy",
        title="Accuracy Curve")

    plot!(p2, val_accs, label="Validation Accuracy")

    display(p1)
    display(p2)
end

# =========================================================
# 🔹 CONFUSION MATRIX
# =========================================================
function confusion_matrix_plot(y_true, y_pred)

    tp = sum((y_pred .== 1) .& (y_true .== 1))
    tn = sum((y_pred .== 0) .& (y_true .== 0))
    fp = sum((y_pred .== 1) .& (y_true .== 0))
    fn = sum((y_pred .== 0) .& (y_true .== 1))

    mat = [tn fp;
           fn tp]

    p = heatmap(mat,
        xticks=(1:2, ["Pred 0", "Pred 1"]),
        yticks=(1:2, ["True 0", "True 1"]),
        title="Confusion Matrix",
        c=:blues,
        annot=true)

    display(p)
end

# =========================================================
# 🔹 F1 SCORE
# =========================================================
function f1_score(y_true, y_pred)

    tp = sum((y_pred .== 1) .& (y_true .== 1))
    fp = sum((y_pred .== 1) .& (y_true .== 0))
    fn = sum((y_pred .== 0) .& (y_true .== 1))

    precision = tp / (tp + fp + 1e-8)
    recall    = tp / (tp + fn + 1e-8)

    f1 = 2 * (precision * recall) / (precision + recall + 1e-8)

    println("Precision: ", precision)
    println("Recall: ", recall)
    println("F1-score: ", f1)

    return f1
end

# =========================================================
# 🔹 ROC CURVE
# =========================================================
function roc_curve_plot(y_true, y_scores)

    thresholds = 0:0.01:1

    tpr_list = Float64[]
    fpr_list = Float64[]

    for t in thresholds

        y_pred = y_scores .>= t

        tp = sum((y_pred .== 1) .& (y_true .== 1))
        fp = sum((y_pred .== 1) .& (y_true .== 0))
        tn = sum((y_pred .== 0) .& (y_true .== 0))
        fn = sum((y_pred .== 0) .& (y_true .== 1))

        tpr = tp / (tp + fn + 1e-8)
        fpr = fp / (fp + tn + 1e-8)

        push!(tpr_list, tpr)
        push!(fpr_list, fpr)
    end

    p = plot(fpr_list, tpr_list,
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        title="ROC Curve",
        label="Model",
        legend=:bottomright)

    plot!(p, [0,1], [0,1], linestyle=:dash, label="Random")

    display(p)
end

# =========================================================
# 🔹 PRECISION-RECALL CURVE
# =========================================================
function precision_recall_plot(y_true, y_scores)

    thresholds = 0:0.01:1

    precisions = Float64[]
    recalls = Float64[]

    for t in thresholds

        y_pred = y_scores .>= t

        tp = sum((y_pred .== 1) .& (y_true .== 1))
        fp = sum((y_pred .== 1) .& (y_true .== 0))
        fn = sum((y_pred .== 0) .& (y_true .== 1))

        precision = tp / (tp + fp + 1e-8)
        recall    = tp / (tp + fn + 1e-8)

        push!(precisions, precision)
        push!(recalls, recall)
    end

    p = plot(recalls, precisions,
        xlabel="Recall",
        ylabel="Precision",
        title="Precision-Recall Curve",
        label="Model",
        legend=:bottomleft)

    display(p)
end

end