module Train

using Flux
using Statistics
using Random
using Optimisers

export train_model, binary_accuracy

function binary_accuracy(ŷ, y)
    preds = ŷ .> 0.5
    return mean(preds .== y)
end

# =========================================================
# AUGMENTATION — uniquement pour mode :custom
# =========================================================

function augment(img::Array{Float32,2})
    if rand() > 0.5
        img = img[:, end:-1:1]
    end
    if rand() > 0.5
        img = img[end:-1:1, :]
    end
    img = img .+ Float32.(randn(size(img)) .* 0.02)
    img = clamp.(img, 0f0, 1f0)
    return img
end

# =========================================================
# MAKE BATCH — gère les deux formats
#
# mode :custom → images sont Array{Float32,2} shape (H, W)
#   batch shape : (W, H, 1, N)
#
# mode :resnet → images sont Array{Float32,3} shape (H, W, 3)
#   batch shape : (W, H, 3, N)
# =========================================================

function make_batch(images, labels, idxs; augment_data=false, model_type=:custom)
    if model_type == :custom
        xs = map(idxs) do i
            img = images[i]   # (H, W)
            if augment_data
                img = augment(img)
            end
            reshape(img, size(img,2), size(img,1), 1)   # → (W, H, 1)
        end
    else  # :resnet
        xs = map(idxs) do i
            img = images[i]   # (H, W, 3)
            # permute vers (W, H, 3) — convention Flux
            permutedims(img, (2, 1, 3))
        end
    end

    xs = cat(xs...; dims=4)
    ys = reshape(Float32.(labels[idxs]), 1, :)
    return xs, ys
end

# =========================================================
# DIAGNOSTIC
# =========================================================

function diagnose_data(train_images, train_labels, model_type)
    println("\n=== DIAGNOSTIC DONNÉES ===")
    println("Nb images train   : ", length(train_images))
    println("Shape image[1]    : ", size(train_images[1]))
    println("Mode              : ", model_type == :custom ? "grayscale (Custom)" : "RGB (ResNet)")
    ratio = sum(train_labels .== 1) / length(train_labels)
    balanced = abs(ratio - 0.5) < 0.1
    println("Ratio pos/total   : ", round(ratio, digits=3), balanced ? " ✓ équilibré" : " ⚠ déséquilibré")
    println("==========================\n")
end

# =========================================================
# ENTRAÎNEMENT
# =========================================================

function train_model(train_data, val_data, model;
                     epochs=30,
                     batch_size=32,
                     lr=1e-3,
                     patience=7,
                     use_augmentation=false,
                     model_type=:custom)

    train_images, train_labels = train_data
    val_images,   val_labels   = val_data

    diagnose_data(train_images, train_labels, model_type)

    if use_augmentation
        println("Data augmentation : ✓ activée")
    else
        println("Data augmentation : ✗ désactivée")
    end
    println()

    N   = length(train_images)
    opt = Optimisers.Adam(lr)
    st  = Optimisers.setup(opt, model)

    train_losses = Float32[]
    val_losses   = Float32[]
    train_accs   = Float32[]
    val_accs     = Float32[]

    best_val_loss    = Inf32
    patience_counter = 0
    lr_current       = Float32(lr)
    lr_counter       = 0
    lr_patience      = 3

    function safe_bce(ŷ, y)
        ŷ_c = clamp.(ŷ, 1f-7, 1f0 - 1f-7)
        return Flux.binarycrossentropy(ŷ_c, y)
    end

    loss_fn(m, x, y) = safe_bce(m(x), y)

    for epoch in 1:epochs

        Flux.trainmode!(model)

        idxs       = shuffle(1:N)
        epoch_loss = 0f0
        epoch_acc  = 0f0
        nb_batches = 0

        for i in 1:batch_size:N
            batch_idxs = idxs[i:min(i+batch_size-1, N)]
            x_batch, y_batch = make_batch(train_images, train_labels, batch_idxs;
                                          augment_data=use_augmentation,
                                          model_type=model_type)

            gs = gradient(m -> loss_fn(m, x_batch, y_batch), model)[1]
            st, model = Optimisers.update(st, model, gs)

            Flux.testmode!(model)
            ŷ = model(x_batch)
            epoch_loss += safe_bce(ŷ, y_batch)
            epoch_acc  += binary_accuracy(ŷ, y_batch)
            Flux.trainmode!(model)

            nb_batches += 1
        end

        train_loss = epoch_loss / nb_batches
        train_acc  = epoch_acc  / nb_batches

        Flux.testmode!(model)

        val_loss = 0f0
        val_acc  = 0f0
        nb_val   = 0

        for i in 1:batch_size:length(val_images)
            batch_idxs = i:min(i+batch_size-1, length(val_images))
            x_batch, y_batch = make_batch(val_images, val_labels, batch_idxs;
                                          augment_data=false,
                                          model_type=model_type)
            ŷ = model(x_batch)
            val_loss += safe_bce(ŷ, y_batch)
            val_acc  += binary_accuracy(ŷ, y_batch)
            nb_val   += 1
        end

        val_loss /= nb_val
        val_acc  /= nb_val

        push!(train_losses, train_loss)
        push!(val_losses,   val_loss)
        push!(train_accs,   train_acc)
        push!(val_accs,     val_acc)

        println("Epoch $epoch | Train Loss: $(round(train_loss, digits=4)) Acc: $(round(train_acc*100, digits=1))% | Val Loss: $(round(val_loss, digits=4)) Acc: $(round(val_acc*100, digits=1))%")

        if val_loss < best_val_loss
            best_val_loss    = val_loss
            patience_counter = 0
            lr_counter       = 0
        else
            patience_counter += 1
            lr_counter       += 1
        end

        if lr_counter >= lr_patience
            lr_current *= 0.5f0
            opt = Optimisers.Adam(lr_current)
            st  = Optimisers.setup(opt, model)
            println("  → LR réduit à $lr_current")
            lr_counter = 0
        end

        if patience_counter >= patience
            println("Early stopping à l'epoch $epoch (best val loss: $(round(best_val_loss, digits=4)))")
            break
        end
    end

    return train_losses, val_losses, train_accs, val_accs
end

end