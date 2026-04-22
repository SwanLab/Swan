module Preprocess

using Images                    #manipuler les images
using ImageTransformations      #transformer les images
using FileIO                    #charger les fichiers
using Random                    #mélanger les données

export load_dataset, split_dataset

# -----------------------------
# 🔹 Charger une image
# -----------------------------
function load_image(path::String, img_size::Tuple{Int, Int})
    img = load(path)

    # Convertir en grayscale
    img_gray = Gray.(img)

    # Resize
    img_resized = imresize(img_gray, img_size)   #passe toutes les images à la meme taille

    # Convertir en Float32 array
    img_array = Float32.(channelview(img_resized))

    return img_array
end

# -----------------------------
# 🔹 Charger dataset complet
# -----------------------------
function load_dataset(root_dir::String; img_size=(128,128), max_per_class=nothing)

    images = []
    labels = Int[]                 # 0 ou 1, pas de fissure ou fissure

    pos_dir = joinpath(root_dir, "Positive")
    neg_dir = joinpath(root_dir, "Negative")

    pos_files = readdir(pos_dir)
    neg_files = readdir(neg_dir)

    if max_per_class !== nothing
        pos_files = pos_files[1:min(end, max_per_class)]
        neg_files = neg_files[1:min(end, max_per_class)]
    end

    # 🔹 Positives (label = 1)
    for f in pos_files
        path = joinpath(pos_dir, f)
        img = load_image(path, img_size)
        push!(images, img)
        push!(labels, 1)
    end

    # 🔹 Négatives (label = 0)
    for f in neg_files
        path = joinpath(neg_dir, f)
        img = load_image(path, img_size)
        push!(images, img)
        push!(labels, 0)
    end

    # Mélange des images P et N
    idx = shuffle(1:length(images))
    images = images[idx]
    labels = labels[idx]

    return images, labels    
end

# -----------------------------
# 🔹 Split train / val / test
# -----------------------------
function split_dataset(images, labels;
                       train_ratio=0.7,
                       val_ratio=0.15)

    N = length(images)

    n_train = Int(floor(train_ratio * N))
    n_val = Int(floor(val_ratio * N))
    n_test = N - n_train - n_val

    train_data = (images[1:n_train], labels[1:n_train])
    val_data = (images[n_train+1:n_train+n_val], labels[n_train+1:n_train+n_val])
    test_data = (images[n_train+n_val+1:end], labels[n_train+n_val+1:end])

    return train_data, val_data, test_data
end

end