# === Chargement des modules ===

using Revise
#include("MyProject.jl")
using .MyProject
using .MyProject.Data
using .MyProject.OptimizationProblemNN
using Plots
using Statistics
using DataFrames
using Random

# === Hyperparamètres ===
pol_deg       = 1
testratio     = 20
λ             = 0.0
learning_rate = 5e-4
hidden_layers = fill(128, 6)

# === Configuration ===
s = Dict{String, Any}()

# Fichier CSV nettoyé (sans header, colonnes numériques uniquement)
# Colonnes : Latitude | Longitude | SpeedOverGround | CourseOverGround |
#            Temperature | Humidity | Pressure | Device1_WindSpeed |
#            Device1_WindAngle | Device1_TrueWindAngle | Device1_TrueWindSpeed |
#            Device1_TrueWindDirection | rpm
s["fileName"] = joinpath(@__DIR__, "Datasets", "merged_clean.csv")

s["polynomialOrder"] = pol_deg
s["testRatio"]       = testratio

s["networkParams"] = Dict(
    "hiddenLayers" => hidden_layers,
    "HUtype"       => "ReLU",
    "OUtype"       => "linear"
)

s["optimizerParams"] = Dict(
    "learningRate" => learning_rate,
    "maxEpochs"    => 3000,
    "type"         => "Adam"
)

s["costParams"] = Dict(
    "λ"        => λ,
    "costType" => "L2"
)

# Inputs  : tout sauf SpeedOverGround (col 3)
# Output  : SpeedOverGround (col 3)
# Col 9 dans X = Device1_TrueWindAngle → transformée en cosd() dans Data.jl
s["xFeatures"] = [1, 2, 4, 8, 9, 10, 11, 12]
s["yFeatures"] = [3]

# === Chargement des données ===
data      = init_data(s)
s["data"] = data

# === Entraînement ===
opt          = init_OptimizationProblemNN(s)
optimizer, θ = solve(opt)
plot_cost(optimizer)

# === Données de test ===
Xtest, Ytest = get_test_data(opt)

# === Réseau entraîné ===
network = get_network(opt)

# === Prédictions (espace normalisé) ===
Ypred = compute_output_values(opt, Xtest, θ)

# === Histogrammes (espace normalisé) ===
edges = range(-3, 3, length=31)
display(histogram(vec(Ypred), bins=edges, title="Distribution de Ypred (normalisé)"))
display(histogram(vec(Ytest), bins=edges, title="Distribution de Ytest (normalisé)"))

# === Visualisation prédictions vs réalité ===
plot_predictions(opt, θ)

# === Dénormalisation ===
Ypred_denorm = vec(Ypred) .* data.sigmaY[1] .+ data.muY[1]
Ytest_denorm = vec(Ytest) .* data.sigmaY[1] .+ data.muY[1]
Xtest_denorm = Xtest .* data.sigmaX .+ data.muX

# === MSE ===
mse = mean((Ypred_denorm .- Ytest_denorm) .^ 2)
rmse = sqrt(mse)
println("MSE  sur les données de test : $(round(mse,  digits=6)) (m/s)²")
println("RMSE sur les données de test : $(round(rmse, digits=4)) m/s")

# === Tableau de résultats ===
difference = Ytest_denorm .- Ypred_denorm

input_data = DataFrame(
    Latitude                = Xtest_denorm[:, 1],
    Longitude               = Xtest_denorm[:, 2],
    CourseOverGround        = Xtest_denorm[:, 3],
    Device1_WindSpeed       = Xtest_denorm[:, 4],
    Device1_WindAngle       = Xtest_denorm[:, 5],
    Device1_TrueWindAngle   = Xtest_denorm[:, 6],
    Device1_TrueWindSpeed   = Xtest_denorm[:, 7],
    Device1_TrueWindDir     = Xtest_denorm[:, 8],
)

output_data = DataFrame(
    Speed_real       = Ytest_denorm,
    Speed_prediction = Ypred_denorm,
    Difference       = difference
)

result_table = hcat(input_data, output_data)
println("Tableau de résultats :")
display(result_table)

# === Permutation Importance ===
feature_names = ["Latitude", "Longitude", "CourseOverGround", "WindSpeed",
                 "WindAngle", "TrueWindAngle", "TrueWindSpeed", "TrueWindDir"]

n_features = size(Xtest, 2)
rmse_base  = sqrt(mean((vec(compute_output_values(opt, Xtest, θ)) .* data.sigmaY[1] .+ data.muY[1]
                        .- Ytest_denorm) .^ 2))

importance = zeros(n_features)
for i in 1:n_features
    Xtest_perm       = copy(Xtest)
    Xtest_perm[:, i] = Xtest_perm[randperm(size(Xtest, 1)), i]
    Yperm            = vec(compute_output_values(opt, Xtest_perm, θ)) .* data.sigmaY[1] .+ data.muY[1]
    rmse_perm        = sqrt(mean((Yperm .- Ytest_denorm) .^ 2))
    importance[i]    = rmse_perm - rmse_base
end

bar_order = sortperm(importance, rev=true)
display(bar(feature_names[bar_order], importance[bar_order];
    title    = "Permutation Importance",
    ylabel   = "Augmentation RMSE (m/s)",
    xlabel   = "Feature",
    legend   = false,
    color    = :steelblue,
    rotation = 30))

    # === Gradient-based Sensitivity ===
grads = vcat([compute_gradient(opt, Xtest[i:i, :], θ) for i in 1:size(Xtest, 1)]...)
grad_mean = vec(mean(abs.(grads), dims=1))

bar_order2 = sortperm(grad_mean, rev=true)
display(bar(feature_names[bar_order2], grad_mean[bar_order2];
    title    = "Gradient-based Sensitivity",
    ylabel   = "∂output/∂input (moyenne |gradient|)",
    xlabel   = "Feature",
    legend   = false,
    color    = :darkorange,
    rotation = 30))


# === One-at-a-time (OAT) ===
n_points   = 50
oat_plots  = []

for i in 1:n_features
    # Plage de variation de la feature i (espace normalisé)
    x_min = minimum(Xtest[:, i])
    x_max = maximum(Xtest[:, i])
    x_range = range(x_min, x_max, length=n_points)

    # Echantillon de base : toutes les features à leur moyenne (= 0 normalisé)
    X_oat = zeros(n_points, n_features)

    # On fait varier uniquement la feature i
    X_oat[:, i] = collect(x_range)

    # Prédiction et dénormalisation
    Y_oat = vec(compute_output_values(opt, X_oat, θ)) .* data.sigmaY[1] .+ data.muY[1]

    # Dénormalisation de l'axe x pour affichage lisible
    x_denorm = collect(x_range) .* data.sigmaX[i] .+ data.muX[i]

    push!(oat_plots, plot(x_denorm, Y_oat;
        title    = feature_names[i],
        xlabel   = feature_names[i],
        ylabel   = "Vitesse (m/s)",
        legend   = false,
        linewidth = 2,
        color    = :green))
end

for i in 1:n_features
    display(oat_plots[i])
end

# === Partial Dependence Plots 2D ===
n_grid = 30  # résolution de la grille

function pdp2d(opt, θ, data, Xtest, i, j, n_grid)
    xi_range = range(minimum(Xtest[:, i]), maximum(Xtest[:, i]), length=n_grid)
    xj_range = range(minimum(Xtest[:, j]), maximum(Xtest[:, j]), length=n_grid)

    Z = zeros(n_grid, n_grid)
    for (ci, xi) in enumerate(xi_range)
        for (cj, xj) in enumerate(xj_range)
            X_base        = zeros(1, size(Xtest, 2))
            X_base[1, i]  = xi
            X_base[1, j]  = xj
            Z[cj, ci]     = compute_output_values(opt, X_base, θ)[1] * data.sigmaY[1] + data.muY[1]
        end
    end

    # Dénormalisation des axes
    xi_denorm = collect(xi_range) .* data.sigmaX[i] .+ data.muX[i]
    xj_denorm = collect(xj_range) .* data.sigmaX[j] .+ data.muX[j]

    return xi_denorm, xj_denorm, Z
end

pairs = [
    (4, 6, "WindSpeed",    "TrueWindAngle"),
    (4, 7, "WindSpeed",    "TrueWindSpeed"),
    (4, 8, "WindSpeed",    "TrueWindDir"),
]

pdp_plots = []
for (i, j, name_i, name_j) in pairs
    xi_d, xj_d, Z = pdp2d(opt, θ, data, Xtest, i, j, n_grid)
    push!(pdp_plots, heatmap(xi_d, xj_d, Z;
        xlabel   = name_i,
        ylabel   = name_j,
        title    = "$name_i × $name_j",
        color    = :viridis,
        colorbar_title = "Vitesse (m/s)"))
end

display(plot(pdp_plots...,
    layout = (1, 3),
    size   = (1400, 450)))