# === Chargement des modules ===

using Revise
#include("MyProject.jl")
using .MyProject
#using .OptimizationProblemNN
#using .Data
using .MyProject.Data
using .MyProject.OptimizationProblemNN
using Plots
using Statistics
using DataFrames

# === Hyperparamètres ===
pol_deg       = 1
testratio     = 20
λ             = 0.0
learning_rate = 1e-3      # learning rate standard pour Adam
hidden_layers = fill(128, 6)

# === Configuration ===
s = Dict{String, Any}()
s["fileName"] = joinpath(@__DIR__, "Datasets", "Resultados2.csv")

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
    "type"         => "Adam" # SGD or Adam
)

s["costParams"] = Dict(
    "λ"        => λ,
    "costType" => "L2"
)

s["xFeatures"] = [1, 2, 3, 5, 6, 7]
s["yFeatures"] = [4]

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
edges = range(-1, 2, length=31)
display(histogram(Ypred, bins=edges, title="Distribution de Ypred (normalisé)"))
display(histogram(Ytest, bins=edges, title="Distribution de Ytest (normalisé)"))

# === Visualisation prédictions vs réalité ===
plot_predictions(opt, θ)

# === Dénormalisation ===
Ypred = Ypred .* data.sigmaY .+ data.muY
Ytest = Ytest .* data.sigmaY .+ data.muY
Xtest = Xtest .* data.sigmaX .+ data.muX

# === MSE ===
mse = mean((Ypred .- Ytest) .^ 2)
println("MSE sur les données de test : $(round(mse, digits=6))")

# === Tableau de résultats ===
difference  = Ytest .- Ypred
input_data  = DataFrame(Xtest, [:rpm, :Windy_cosine, :Windy_ms, :Speed3, :Yaw, :Pitch, :Roll])
output_data = DataFrame(
    Cons_real       = vec(Ytest),
    Cons_prediction = vec(Ypred),
    Difference      = vec(difference)
)

result_table = hcat(input_data, output_data)
println("Tableau de résultats :")
display(result_table)