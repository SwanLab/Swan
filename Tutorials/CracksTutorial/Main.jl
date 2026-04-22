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
# 🔹 PARAMÈTRES GLOBAUX
# =========================================================

root_dir = "C:/Users/couwa/Documents/Stage_CIMNE/Classification/archive"  # <-- à modifier

img_size = (128, 128)

batch_size = 32
epochs = 20
learning_rate = 1e-3
max_images = 200   # pour test rapide

# =========================================================
# 🔹 DATA
# =========================================================

println("Loading dataset...")

images, labels = load_dataset(root_dir;
                              img_size=img_size,
                              max_per_class=max_images)

train_data, val_data, test_data = split_dataset(images, labels)

println("Dataset loaded:")
println("Train: ", length(train_data[1]))
println("Val:   ", length(val_data[1]))
println("Test:  ", length(test_data[1]))

# =========================================================
# 🔹 MODEL
# =========================================================

println("Building model...")

model = build_model()

# =========================================================
# 🔹 TRAINING
# =========================================================

println("Training started...")

train_losses, val_losses, train_accs, val_accs =
    train_model(model, train_data, val_data;
                epochs=epochs,
                batch_size=batch_size,
                lr=learning_rate,
                patience=5)

# =========================================================
# 🔹 PLOTS
# =========================================================

println("Plotting results...")

plot_metrics(train_losses, val_losses,
             train_accs, val_accs)

# =========================================================
# 🔹 TEST FINAL
# =========================================================

println("Final evaluation on test set...")

test_images, test_labels = test_data

# predictions proba
y_scores = Float32[]

for x in test_images
    x = reshape(Float32.(x), size(x,1), size(x,2), 1, 1)
    ŷ = model(x)
    push!(y_scores, ŷ[1])   # scalaire propre
end

println("Scores sample: ", y_scores[1:10])
println("Mean score: ", mean(y_scores))
println("Min score: ", minimum(y_scores))
println("Max score: ", maximum(y_scores))

# binary predictions
y_pred = y_scores .> 0.5

# metrics avancées
f1_score(test_labels, y_pred)
confusion_matrix_plot(test_labels, y_pred)
roc_curve_plot(test_labels, y_scores)
precision_recall_plot(test_labels, y_scores)

println("Done.")