module Preprocess

using Images
using ImageTransformations
using FileIO
using Random

export load_dataset, load_dataset_resnet, split_dataset

# =========================================================
# MODE CUSTOM — grayscale, normalisation min-max
# =========================================================

function load_image_custom(path::String, img_size::Tuple{Int,Int})
    img      = load(path)
    img_gray = Gray.(img)
    img_res  = imresize(img_gray, img_size)
    img_arr  = Float32.(channelview(img_res))

    # Normalisation min-max par image : étire dans [0,1]
    lo, hi  = minimum(img_arr), maximum(img_arr)
    img_arr = (img_arr .- lo) ./ (hi - lo + 1f-8)

    return img_arr   # shape (H, W)
end

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
        push!(images, load_image_custom(path, img_size))
        push!(labels, 1)
    end

    println("Loading $(length(neg_files)) negative images...")
    for path in neg_files
        push!(images, load_image_custom(path, img_size))
        push!(labels, 0)
    end

    idx    = shuffle(1:length(images))
    images = images[idx]
    labels = labels[idx]

    println("Dataset: $(length(images)) images ($(sum(labels)) positives, $(sum(labels.==0)) negatives)")
    return images, labels
end

# =========================================================
# MODE RESNET — RGB, normalisation ImageNet, 128×128 fixe
#
# ResNet18 a été entraîné sur ImageNet avec ces stats :
#   mean = [0.485, 0.456, 0.406]  par channel R, G, B
#   std  = [0.229, 0.224, 0.225]  par channel R, G, B
#
# Appliquer la même normalisation permet aux activations
# internes de ResNet de rester dans leur plage habituelle,
# ce qui rend le transfer learning beaucoup plus efficace.
#
# Format de sortie : Array{Float32, 3} de shape (H, W, 3)
# =========================================================

const RESNET_IMG_SIZE = (128, 128)

# Moyennes et écarts-types ImageNet par channel R, G, B
const IMAGENET_MEAN = Float32[0.485, 0.456, 0.406]
const IMAGENET_STD  = Float32[0.229, 0.224, 0.225]

function load_image_resnet(path::String)
    img     = load(path)
    img_res = imresize(img, RESNET_IMG_SIZE)

    # Convertit en RGB pour garantir 3 channels
    # (certaines images peuvent être RGBA ou grayscale)
    img_rgb = RGB.(img_res)

    # Extrait les 3 channels : shape (3, H, W)
    img_arr = Float32.(channelview(img_rgb))
    # img_arr[1, :, :] = channel R
    # img_arr[2, :, :] = channel G
    # img_arr[3, :, :] = channel B

    # Normalisation ImageNet : (pixel - mean) / std par channel
    for c in 1:3
        img_arr[c, :, :] = (img_arr[c, :, :] .- IMAGENET_MEAN[c]) ./ IMAGENET_STD[c]
    end

    # Retourne (H, W, 3) — convention plus lisible, make_batch gère la suite
    return permutedims(img_arr, (2, 3, 1))
end

function load_dataset_resnet(root_dir::String; max_per_class=nothing)

    images = Vector{Array{Float32,3}}()
    labels = Int[]

    pos_dir = joinpath(root_dir, "Positive")
    neg_dir = joinpath(root_dir, "Negative")

    pos_files = readdir(pos_dir, join=true)
    neg_files = readdir(neg_dir, join=true)

    if max_per_class !== nothing
        pos_files = pos_files[1:min(length(pos_files), max_per_class)]
        neg_files = neg_files[1:min(length(neg_files), max_per_class)]
    end

    println("Loading $(length(pos_files)) positive images (ResNet / RGB 128×128)...")
    for path in pos_files
        push!(images, load_image_resnet(path))
        push!(labels, 1)
    end

    println("Loading $(length(neg_files)) negative images (ResNet / RGB 128×128)...")
    for path in neg_files
        push!(images, load_image_resnet(path))
        push!(labels, 0)
    end

    idx    = shuffle(1:length(images))
    images = images[idx]
    labels = labels[idx]

    println("Dataset: $(length(images)) images ($(sum(labels)) positives, $(sum(labels.==0)) negatives)")
    return images, labels
end

# =========================================================
# SPLIT STRATIFIÉ PAR BLOCS — commun aux deux modes
#
# Les images consécutives dans ce dataset viennent du même
# mur physique → on découpe en blocs contigus pour éviter
# que des images quasi-identiques se retrouvent dans train
# et test simultanément (data leakage).
# =========================================================

function split_dataset(images, labels;
                       train_ratio=0.7,
                       val_ratio=0.15,
                       block_size=50)

    pos_idx = findall(==(1), labels)
    neg_idx = findall(==(0), labels)

    function block_split(idx, r_train, r_val)
        n      = length(idx)
        blocks = [idx[i:min(i+block_size-1, n)] for i in 1:block_size:n]
        shuffle!(blocks)

        n_blocks       = length(blocks)
        n_train_blocks = Int(floor(r_train * n_blocks))
        n_val_blocks   = Int(floor(r_val   * n_blocks))

        train_idx = vcat(blocks[1:n_train_blocks]...)
        val_idx   = vcat(blocks[n_train_blocks+1:n_train_blocks+n_val_blocks]...)
        test_idx  = vcat(blocks[n_train_blocks+n_val_blocks+1:end]...)

        return train_idx, val_idx, test_idx
    end

    pos_tr, pos_v, pos_te = block_split(pos_idx, train_ratio, val_ratio)
    neg_tr, neg_v, neg_te = block_split(neg_idx, train_ratio, val_ratio)

    train_idx = shuffle(vcat(pos_tr, neg_tr))
    val_idx   = shuffle(vcat(pos_v,  neg_v))
    test_idx  = shuffle(vcat(pos_te, neg_te))

    train_data = (images[train_idx], labels[train_idx])
    val_data   = (images[val_idx],   labels[val_idx])
    test_data  = (images[test_idx],  labels[test_idx])

    println("Split par blocs de $block_size | Train: $(length(train_idx)) | Val: $(length(val_idx)) | Test: $(length(test_idx))")
    return train_data, val_data, test_data
end

end