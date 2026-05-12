# =========================================================
# visualize_preprocessing.jl
#
# Script standalone — à lancer séparément du pipeline principal
# But : afficher les images TELLES QUE LE CNN LES REÇOIT
#       après grayscale + resize + normalisation min-max
#       et optionnellement avec la data augmentation
#
# Usage :
#   julia visualize_preprocessing.jl
# =========================================================

using Images
using ImageTransformations
using FileIO
using Plots
using Statistics
gr()

# =========================================================
# PARAMÈTRES
# =========================================================

root_dir = joinpath(@__DIR__, "archive")
img_size = (128, 128)   # doit correspondre à Main.jl

n_samples = 8           # total affiché (moitié pos, moitié neg)

output_path = joinpath(@__DIR__, "plots", "preprocessed_samples.png")

# --- Options d'augmentation ---
# Chaque option est indépendante : active ce que tu veux visualiser
use_flip_h       = false   # flip horizontal
use_flip_v       = false   # flip vertical
use_gaussian_noise = true # bruit gaussien (simule texture béton)
use_gaussian_blur  = false # flou gaussien (simule mise au point imparfaite)

gaussian_noise_std = 0.02  # écart-type du bruit  (0.01=léger, 0.05=fort)
gaussian_blur_std  = 1.0   # écart-type du flou   (0.5=léger, 2.0=fort)

# =========================================================
# Preprocessing — reproduction exacte de module_Preprocess
# =========================================================

function load_and_preprocess(path::String, img_size::Tuple{Int,Int})
    img      = load(path)
    img_gray = Gray.(img)
    img_res  = imresize(img_gray, img_size)
    img_arr  = Float32.(channelview(img_res))

    # Normalisation min-max
    lo, hi  = minimum(img_arr), maximum(img_arr)
    img_arr = (img_arr .- lo) ./ (hi - lo + 1f-8)

    return img_arr
end

# =========================================================
# Augmentation — même logique que module_Train
# =========================================================

# Noyau gaussien 1D pour le flou
function gaussian_kernel_1d(size::Int, sigma::Float64)
    center = size ÷ 2 + 1
    k = [exp(-((i - center)^2) / (2 * sigma^2)) for i in 1:size]
    return Float32.(k ./ sum(k))
end

# Convolution 1D appliquée sur chaque ligne puis chaque colonne (flou séparable)
function apply_gaussian_blur(img::Array{Float32,2}, sigma::Float64)
    ksize  = max(3, 2 * Int(ceil(3 * sigma)) + 1)   # taille du noyau : 3σ de chaque côté
    kernel = gaussian_kernel_1d(ksize, sigma)
    pad    = ksize ÷ 2

    H, W = size(img)
    out  = copy(img)

    # Flou horizontal (sur chaque ligne)
    for i in 1:H
        for j in 1:W
            val = 0f0
            for k in 1:ksize
                jj = clamp(j - pad + k - 1, 1, W)
                val += kernel[k] * img[i, jj]
            end
            out[i, j] = val
        end
    end

    tmp = copy(out)

    # Flou vertical (sur chaque colonne)
    for i in 1:H
        for j in 1:W
            val = 0f0
            for k in 1:ksize
                ii = clamp(i - pad + k - 1, 1, H)
                val += kernel[k] * tmp[ii, j]
            end
            out[i, j] = val
        end
    end

    return out
end

function augment(img::Array{Float32,2};
                 flip_h=false, flip_v=false,
                 add_noise=false, noise_std=0.02,
                 add_blur=false,  blur_std=1.0)

    if flip_h
        img = img[:, end:-1:1]
    end

    if flip_v
        img = img[end:-1:1, :]
    end

    if add_noise
        img = img .+ Float32.(randn(size(img)) .* noise_std)
        img = clamp.(img, 0f0, 1f0)
    end

    if add_blur
        img = apply_gaussian_blur(img, blur_std)
    end

    return img
end

# =========================================================
# Chargement
# =========================================================

pos_dir   = joinpath(root_dir, "Positive")
neg_dir   = joinpath(root_dir, "Negative")

pos_files = readdir(pos_dir, join=true)[1:n_samples÷2]
neg_files = readdir(neg_dir, join=true)[1:n_samples÷2]

# Résumé des options actives
aug_active = any([use_flip_h, use_flip_v, use_gaussian_noise, use_gaussian_blur])
println("\n=== Options d'augmentation ===")
println("Flip horizontal  : ", use_flip_h       ? "✓" : "✗")
println("Flip vertical    : ", use_flip_v       ? "✓" : "✗")
println("Bruit gaussien   : ", use_gaussian_noise ? "✓ (std=$(gaussian_noise_std))" : "✗")
println("Flou gaussien    : ", use_gaussian_blur  ? "✓ (std=$(gaussian_blur_std))"  : "✗")
println()

# =========================================================
# Construction des plots
# =========================================================

plots_list = []

for (label, files, tag) in [("POS", pos_files, "fissure"), ("NEG", neg_files, "sain")]
    for (i, path) in enumerate(files)
        img_arr = load_and_preprocess(path, img_size)

        if aug_active
            img_arr = augment(img_arr;
                flip_h     = use_flip_h,
                flip_v     = use_flip_v,
                add_noise  = use_gaussian_noise,
                noise_std  = gaussian_noise_std,
                add_blur   = use_gaussian_blur,
                blur_std   = gaussian_blur_std)
        end

        println("$label $i — Min: $(round(minimum(img_arr), digits=3))  Max: $(round(maximum(img_arr), digits=3))  Mean: $(round(mean(img_arr), digits=3))")

        p = heatmap(img_arr,
            color=:grays,
            axis=false,
            colorbar=false,
            title="$label $i",
            titlefontsize=8,
            aspect_ratio=:equal)

        push!(plots_list, p)
    end
end

# Titre dynamique selon les options actives
aug_labels = String[]
use_flip_h         && push!(aug_labels, "flip H")
use_flip_v         && push!(aug_labels, "flip V")
use_gaussian_noise && push!(aug_labels, "bruit σ=$(gaussian_noise_std)")
use_gaussian_blur  && push!(aug_labels, "flou σ=$(gaussian_blur_std)")

title_aug = isempty(aug_labels) ? "sans augmentation" : join(aug_labels, " + ")
plot_title = "Images après preprocessing — $(title_aug)\n(ligne 1 : fissures | ligne 2 : sans fissure)"

n_cols = n_samples ÷ 2
fig = plot(plots_list...,
    layout=(2, n_cols),
    size=(n_cols * 200, 500),
    plot_title=plot_title,
    plot_titlefontsize=9)

isdir(dirname(output_path)) || mkdir(dirname(output_path))
savefig(fig, output_path)
println("\n→ Saved: $output_path")