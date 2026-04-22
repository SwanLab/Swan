module Train

using Flux
using Statistics
using Random
using Optimisers

export train_model, binary_accuracy

# -----------------------------
# 🔹 Accuracy (batch)
# -----------------------------
function binary_accuracy(ŷ, y)
    preds = ŷ .> 0.5
    return mean((preds .== y))
end

# -----------------------------
# 🔹 Préparer un batch
# -----------------------------
function make_batch(images, labels, idxs)

    xs = [Float32.(images[i]) for i in idxs]

    # reshape en (H, W, C)
    xs = [reshape(x, size(x,1), size(x,2), 1) for x in xs]

    # batch = 4e dimension
    xs = cat(xs...; dims=4)

    ys = Float32.(labels[idxs])
    ys = reshape(ys, 1, :)

    return xs, ys
end

# -----------------------------
# 🔹 Entraînement
# -----------------------------
function train_model(model, train_data, val_data;
                     epochs=20,
                     batch_size=32,
                     lr=1e-3,
                     patience=5)

    opt = Optimisers.Adam(lr)
    st = Optimisers.setup(opt, model)

    train_images, train_labels = train_data
    val_images, val_labels = val_data

    N = length(train_images)

    # 🔹 Historique
    train_losses = Float32[]
    val_losses = Float32[]
    train_accs = Float32[]
    val_accs = Float32[]

    # 🔹 Early stopping
    best_val_loss = Inf
    patience_counter = 0

    # -----------------------------
    # 🔹 LOSS FUNCTION (NEW STYLE)
    # -----------------------------
    loss_fn(m, x, y) = Flux.binarycrossentropy(m(x), y)

    for epoch in 1:epochs

        idxs = shuffle(1:N)

        epoch_loss = 0.0
        epoch_acc = 0.0
        nb_batches = 0

        # -----------------------------
        # 🔹 TRAIN
        # -----------------------------
        for i in 1:batch_size:N
            batch_idxs = idxs[i:min(i+batch_size-1, N)]

            x_batch, y_batch = make_batch(train_images, train_labels, batch_idxs)

            loss(m) = loss_fn(m, x_batch, y_batch)

            gs = gradient(loss, model)[1]

            st, model = Optimisers.update(st, model, gs)

            # Metrics
            ŷ = model(x_batch)
            epoch_loss += Flux.binarycrossentropy(ŷ, y_batch)
            epoch_acc += binary_accuracy(ŷ, y_batch)

            nb_batches += 1
        end

        train_loss = epoch_loss / nb_batches
        train_acc = epoch_acc / nb_batches

        # -----------------------------
        # 🔹 VALIDATION
        # -----------------------------
        val_loss = 0.0
        val_acc = 0.0
        nb_val_batches = 0

        for i in 1:batch_size:length(val_images)
            batch_idxs = i:min(i+batch_size-1, length(val_images))

            x_batch, y_batch = make_batch(val_images, val_labels, batch_idxs)

            ŷ = model(x_batch)

            val_loss += Flux.binarycrossentropy(ŷ, y_batch)
            val_acc += binary_accuracy(ŷ, y_batch)

            nb_val_batches += 1
        end

        val_loss /= nb_val_batches
        val_acc /= nb_val_batches

        # -----------------------------
        # 🔹 Stockage
        # -----------------------------
        push!(train_losses, train_loss)
        push!(val_losses, val_loss)
        push!(train_accs, train_acc)
        push!(val_accs, val_acc)

        println("Epoch $epoch")
        println("Train Loss: $train_loss | Train Acc: $train_acc")
        println("Val Loss:   $val_loss   | Val Acc:   $val_acc")
        println("--------------------------------------------------")

        # -----------------------------
        # 🔥 Early Stopping
        # -----------------------------
        if val_loss < best_val_loss
            best_val_loss = val_loss
            patience_counter = 0
        else
            patience_counter += 1
        end

        if patience_counter >= patience
            println(" Early stopping triggered at epoch $epoch")
            break
        end
    end

    return train_losses, val_losses, train_accs, val_accs
end

end