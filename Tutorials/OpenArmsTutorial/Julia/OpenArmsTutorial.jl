# === Include all modules in dependency order ===

include("LearnableVariables/LearnableVariables.jl")  # No dependencies
include("Data/Data.jl")                              # Uses LearnableVariables indirectly (via network)
include("Network/Network.jl")                        # Uses LearnableVariables
include("LossFunctional/LossFunctional.jl")          # Uses Network, Data, LearnableVariables
include("Sh_Func_L2norm/Sh_Func_L2norm.jl")          # Uses LearnableVariables
include("CostNN/CostNN.jl")                          # Uses LossFunctional, Sh_Func_L2norm
include("PlotterNN/PlotterNN.jl")                    # Uses Network, CostNN, Data

# Trainer and optimization methods
include("Trainer/Trainer.jl")                        # Depends on nothing
include("Trainer/SGD/SGD.jl")                                # Uses Trainer
include("Trainer/Nesterov/Nesterov.jl")                      # Uses Trainer
include("Trainer/RMSProp/RMSProp.jl")                        # Uses Trainer
include("Trainer/Fminunc/Fminunc.jl")                        # Uses Trainer

# Final controller
include("OptimizationProblemNN/OptimizationProblemNN.jl")  # Uses everything above


using .OptimizationProblemNN
using .Data
using Plots
using Statistics
using DataFrames
# Tools to see where time is spent and improve performance
using Profile 
using ProfileView

# Initialization of hyperparameters
pol_deg       = 1
testratio     = 20
lambda        = 0.0
learningRate  = 0.2
hiddenLayers  = fill(128, 6)

# INITIALIZATION
# Store dataset file name
s = Dict{String, Any}()
#s["fileName"] = "Resultados2.csv"
s["fileName"] = joinpath(@__DIR__, "Datasets", "Resultados2.csv")

# Load model parameters
s["polynomialOrder"] = pol_deg
s["testRatio"] = testratio
s["networkParams"] = Dict(
    "hiddenLayers" => hiddenLayers,
    "HUtype" => "ReLU",
    "OUtype" => "linear"
)

s["optimizerParams"] = Dict(
    "learningRate" => learningRate,
    "maxEpochs" => 10 # adjust to 10 for fast runs if needed
)

s["costParams"] = Dict(
    "lambda" => lambda,
    "costType" => "L2"
)

# Select the model's features
s["xFeatures"] = [1,2,3,4,5,6,7]
s["yFeatures"] = [8]

# Load data
# data = cHomogData(s)
# data = JuliaData(s)
data = Data.DataStruct(s)
s["data"] = data

# Train the model
opt = OptimizationProblemNN.OptimizationProblemNNStruct(s)
OptimizationProblemNN.solve(opt)
OptimizationProblemNN.plotCostFnc(opt)


# Obtain test data
Xtest, Ytest = OptimizationProblemNN.getTestData(opt)

# Load the trained neural network
network = OptimizationProblemNN.getNetwork(opt)

# Vectorized prediction of Ytest
Ypred = OptimizationProblemNN.computeOutputValues(opt, Xtest)

# Histogram for the distribution of Ypred
edges = range(-1, 2, length=31)  # 30 bins between -1 and 2
#=
hist1 = histogram(Ypred, bins=edges, title="Distribution of predicted Ytest")
display(hist1)

# Histogram for the distribution of Ytest
hist2 = histogram(Ytest, bins=edges, title="Distribution of Test Y")
display(hist2)
=#
# Denormalization
Xtest = Xtest .* data.sigmaX .+ data.muX
Ypred = Ypred .* data.sigmaY .+ data.muY
Ytest = Ytest .* data.sigmaY .+ data.muY

mse = mean((Ypred .- Ytest).^2)
println("Error cuadrático medio (MSE) en los datos de prueba: $(round(mse, digits=6))")

# Compute difference between real and predicted values
difference = Ytest .- Ypred

# Convert Xtest to DataFrame and rename columns
input_data = DataFrame(Xtest, [:rpm, :Windy_cosine, :Windy_ms, :Speed3, :Yaw, :Pitch, :Roll])

# Create output DataFrame
output_data = DataFrame(
    Cons_real = vec(Ytest),
    Cons_prediction = vec(Ypred),
    Difference = vec(difference)
)
# Combine input and output into one result table
result_table = hcat(input_data, output_data)

# Display table
println("Tabla de resultados:")
display(result_table)
