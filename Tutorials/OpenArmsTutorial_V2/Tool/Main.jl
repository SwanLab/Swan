# === Chargement des modules ===

using Revise
#include("MyProject.jl")
using .MyProject
using .MyProject.Data
using .MyProject.OptimizationProblemNN
using Plots
using Statistics
using DataFrames

# === Hyperparamètres ===
pol_deg       = 1
testratio     = 20
λ             = 0.0
learning_rate = 1e-3
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
    "maxEpochs"    => 1000,
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