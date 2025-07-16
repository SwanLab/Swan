# === Include all modules through MyProject ===

using Revise
include("MyProject.jl")
using .MyProject
using .OptimizationProblemNN
using .Data
using Plots
using Statistics
using DataFrames

# === Initialization of hyperparameters ===
pol_deg       = 1
testratio     = 20
λ             = 0.0     
learning_rate = 0.2
hidden_layers = fill(128, 6)

# === INITIALIZATION ===

# Store dataset file name
s = Dict{String, Any}()
s["fileName"] = joinpath(@__DIR__, "Datasets", "Resultados2.csv")

# Load model parameters
s["polynomialOrder"] = pol_deg
s["testRatio"] = testratio
s["networkParams"] = Dict(
    "hiddenLayers" => hidden_layers,
    "HUtype" => "ReLU",
    "OUtype" => "linear"
)

s["optimizerParams"] = Dict(
    "learningRate" => learning_rate,
    "maxEpochs" => 1000,  # adjust to 10 for fast runs if needed
    "type" => "SGD"      # specify your optimizer here
)

s["costParams"] = Dict(
    "λ" => λ,
    "costType" => "L2"
)

# Select the model's features
s["xFeatures"] = [1,2,3,4,5,6,7]
s["yFeatures"] = [8]

# === Load data ===
data = Data.init_data(s)
s["data"] = data

# === Train the model ===
opt = init_OptimizationProblemNN(s)
optimizer, θ = solve(opt)
plot_cost(optimizer)

# === Obtain test data ===
Xtest, Ytest = get_test_data(opt)

# === Load the trained neural network ===
network = get_network(opt)

# === Vectorized prediction of Ytest ===
Ypred = compute_output_values(opt, Xtest, θ)

# === Histogram for predicted Ytest distribution ===
edges = range(-1, 2, length=31)  # 30 bins between -1 and 2
hist1 = histogram(Ypred, bins=edges, title="Distribution of predicted Ytest")
display(hist1)

# === Histogram for actual Ytest distribution ===
hist2 = histogram(Ytest, bins=edges, title="Distribution of Test Y")
display(hist2)

# === Denormalize ===
Xtest = Xtest .* data.sigmaX .+ data.muX
Ypred = Ypred .* data.sigmaY .+ data.muY
Ytest = Ytest .* data.sigmaY .+ data.muY

# === Compute and print MSE ===
mse = mean((Ypred .- Ytest).^2)
println("Error cuadrático medio (MSE) en los datos de prueba: $(round(mse, digits=6))")

# === Compute difference between real and predicted values ===
difference = Ytest .- Ypred

# === Convert Xtest to DataFrame and rename columns ===
input_data = DataFrame(Xtest, [:rpm, :Windy_cosine, :Windy_ms, :Speed3, :Yaw, :Pitch, :Roll])

# === Create output DataFrame ===
output_data = DataFrame(
    Cons_real = vec(Ytest),
    Cons_prediction = vec(Ypred),
    Difference = vec(difference)
)

# === Combine input and output into one result table ===
result_table = hcat(input_data, output_data)

# === Display table ===
println("Tabla de resultados:")
display(result_table)