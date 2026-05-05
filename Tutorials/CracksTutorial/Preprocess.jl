module Preprocess

using Images
using ImageTransformations
using FileIO
using Random

export load_dataset, split_dataset

# -----------------------------
# Charger et normaliser une image
#
# PROBLÈME OBSERVÉ : Min pixel = 0.52, Max pixel = 0.96
# Les images sources sont encodées dans un sous-intervalle
# de [0,1] — le modèle ne voit jamais de pixels sombres
# → il ne peut pas apprendre à distinguer les fissures
# (zones sombres) du béton (zones claires).
#
# FIX : normalisation min-max par image
# Chaque image est étirée pour utiliser tout [0,1].
# Un mur sans fissure aura des valeurs proches et uniformes
# → après normalisation, contraste amplifié → le CNN voit mieux.
# -----------------------------
function load_image(path::String, img_size::Tuple{Int, Int})
    img = load(path)

    img_gray    = Gray.(img)
    img_resized = imresize(img_gray, img_size)
    img_array   = Float32.(channelview(img_resized))    #for type compatibility with Flux

    # Normalisation min-max par image : étire les valeurs dans [0,1]
    lo, hi = minimum(img_array), maximum(img_array)
    img_array = (img_array .- lo) ./ (hi - lo + 1f-8)

    return img_array
end

# -----------------------------
# Charger le dataset complet
# root_dir doit contenir "Positive/" et "Negative/"
# -----------------------------
function load_dataset(root_dir::String; img_size=(64,64), max_per_class=nothing)

    images = Vector{Array{Float32,2}}()
    labels = Int[]

    pos_dir = joinpath(root_dir, "Positive")
    neg_dir = joinpath(root_dir, "Negative")

    pos_files = readdir(pos_dir, join=true)
    neg_files = readdir(neg_dir, join=true)

    if max_per_class !== nothing
        pos_files = pos_files[1:min(length(pos_files), max_per_class)]
        neg_files = neg_files[1:min(length(neg_files), max_per_class)]
    end

    println("Loading $(length(pos_files)) positive images...")
    for path in pos_files
        push!(images, load_image(path, img_size))
        push!(labels, 1)
    end

    println("Loading $(length(neg_files)) negative images...")
    for path in neg_files
        push!(images, load_image(path, img_size))
        push!(labels, 0)
    end

    idx    = shuffle(1:length(images))
    images = images[idx]
    labels = labels[idx]

    println("Dataset: $(length(images)) images ($(sum(labels)) positives, $(sum(labels.==0)) negatives)")
    return images, labels
end

# -----------------------------
# Split stratifié — garantit 50/50 pos/neg dans chaque split
#
# Sans stratification, un shuffle maladroit peut mettre
# 80% de positifs dans le val set → val_loss qui mesure
# autre chose que les vraies performances.
# -----------------------------
function split_dataset(images, labels; train_ratio=0.7, val_ratio=0.15)

    pos_idx = findall(==(1), labels)
    neg_idx = findall(==(0), labels)

    function split_idx(idx, r_train, r_val)
        n       = length(idx)
        n_train = Int(floor(r_train * n))
        n_val   = Int(floor(r_val   * n))
        return idx[1:n_train], idx[n_train+1:n_train+n_val], idx[n_train+n_val+1:end]
    end

    pos_tr, pos_v, pos_te = split_idx(pos_idx, train_ratio, val_ratio)
    neg_tr, neg_v, neg_te = split_idx(neg_idx, train_ratio, val_ratio)

    train_idx = shuffle(vcat(pos_tr, neg_tr))
    val_idx   = shuffle(vcat(pos_v,  neg_v))
    test_idx  = shuffle(vcat(pos_te, neg_te))

    train_data = (images[train_idx], labels[train_idx])
    val_data   = (images[val_idx],   labels[val_idx])
    test_data  = (images[test_idx],  labels[test_idx])

    println("Train: $(length(train_idx)) | Val: $(length(val_idx)) | Test: $(length(test_idx))")
    return train_data, val_data, test_data
end

end