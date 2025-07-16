# === Include all modules in dependency order ===

include("LearnableVariables.jl")  # No dependencies
include("Data.jl")                              # Uses LearnableVariables indirectly (via network)
include("Network.jl")                        # Uses LearnableVariables
include("LossFunctional.jl")          # Uses Network, Data, LearnableVariables
include("Sh_Func_L2norm.jl")          # Uses LearnableVariables
include("CostNN.jl")                          # Uses LossFunctional, Sh_Func_L2norm
include("PlotterNN.jl")                    # Uses Network, CostNN, Data

# Trainer and optimization methods
include("Trainer.jl")                        # Depends on nothing
include("SGD.jl")                                # Uses Trainer
include("Nesterov.jl")                      # Uses Trainer
include("RMSProp.jl")                        # Uses Trainer
include("Fminunc.jl")                        # Uses Trainer

# Final controller
include("OptimizationProblemNN.jl")  # Uses everything above

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
    "maxEpochs" => 100,  # adjust to 10 for fast runs if needed
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
θ = solve(opt)
plot_cost(opt)

# === Obtain test data ===
Xtest, Ytest = get_test_data(opt)

# === Load the trained neural network ===
network = get_network(opt)



