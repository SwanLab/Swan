# =========================================================
# visualize_preprocessing.jl
#
# Script standalone — à lancer séparément du pipeline principal
# But : afficher les images TELLES QUE LE CNN LES REÇOIT
#       après grayscale + resize + normalisation min-max
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
# PARAMÈTRES — adapter selon ta config
# =========================================================

root_dir = joinpath(@__DIR__, "archive")
img_size = (128, 128)   # doit correspondre à ce que tu utilises dans Main.jl

n_samples = 8           # nombre d'images à afficher (4 positives + 4 négatives)

output_path = joinpath(@__DIR__, "plots", "preprocessed_samples.png")

# =========================================================
# Reproduction exacte du preprocessing de module_Preprocess
# (copié ici pour ne pas dépendre du module)
# =========================================================

function load_and_preprocess(path::String, img_size::Tuple{Int,Int})
    img       = load(path)
    img_gray  = Gray.(img)
    img_res   = imresize(img_gray, img_size)
    img_arr   = Float32.(channelview(img_res))

    # Normalisation min-max — identique à module_Preprocess
    lo, hi  = minimum(img_arr), maximum(img_arr)
    img_arr = (img_arr .- lo) ./ (hi - lo + 1f-8)

    return img_arr
end

# =========================================================
# Chargement de quelques images positives et négatives
# =========================================================

pos_dir   = joinpath(root_dir, "Positive")
neg_dir   = joinpath(root_dir, "Negative")

pos_files = readdir(pos_dir, join=true)[1:n_samples÷2]
neg_files = readdir(neg_dir, join=true)[1:n_samples÷2]

println("Chargement de $(length(pos_files)) positives et $(length(neg_files)) négatives...")

# =========================================================
# Affichage
#
# Chaque image est affichée en grayscale dans [0,1]
# C'est EXACTEMENT ce que voit la première couche Conv du CNN
# La normalisation min-max amplifie le contraste :
#   - pixel le plus sombre de l'image → 0.0 (noir)
#   - pixel le plus clair             → 1.0 (blanc)
# =========================================================

plots_list = []

for (i, path) in enumerate(pos_files)
    img_arr = load_and_preprocess(path, img_size)

    println("Positive $i — Min: $(round(minimum(img_arr), digits=3))  Max: $(round(maximum(img_arr), digits=3))  Mean: $(round(mean(img_arr), digits=3))")

    p = heatmap(img_arr,
        color=:grays,
        axis=false,
        colorbar=false,
        title="POS $i",
        titlefontsize=8,
        aspect_ratio=:equal)

    push!(plots_list, p)
end

for (i, path) in enumerate(neg_files)
    img_arr = load_and_preprocess(path, img_size)

    println("Negative $i — Min: $(round(minimum(img_arr), digits=3))  Max: $(round(maximum(img_arr), digits=3))  Mean: $(round(mean(img_arr), digits=3))")

    p = heatmap(img_arr,
        color=:grays,
        axis=false,
        colorbar=false,
        title="NEG $i",
        titlefontsize=8,
        aspect_ratio=:equal)

    push!(plots_list, p)
end

# Grille : 2 lignes × (n_samples/2) colonnes
# Ligne 1 : positives (avec fissures)
# Ligne 2 : négatives (sans fissures)
n_cols = n_samples ÷ 2
fig = plot(plots_list...,
    layout=(2, n_cols),
    size=(n_cols * 200, 500),
    plot_title="Images après preprocessing — telles que vues par le CNN\n(ligne 1 : fissures | ligne 2 : sans fissure)",
    plot_titlefontsize=10)

isdir(dirname(output_path)) || mkdir(dirname(output_path))
savefig(fig, output_path)
println("\n→ Saved: $output_path")
println("Ouvre ce fichier pour voir si les fissures sont visibles après preprocessing.")