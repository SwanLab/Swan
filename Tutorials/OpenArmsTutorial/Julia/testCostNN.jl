include("Network/Network.jl")
include("Sh_Func_L2norm/Sh_Func_L2norm.jl")
include("LossFunctional/LossFunctional.jl")
include("CostNN/CostNN.jl")

using LinearAlgebra
using Random
using .Network
using .Sh_Func_L2norm
using .LossFunctional
using .CostNN
using .Network.LearnableVariables

# Dummy "network" and "data" for LossFunc
n_features = 5
n_outputs = 3
n_samples = 20

# Random input matrix: nSamples x nFeatures
Xtrain = randn(n_samples, n_features)
Ytrain = zeros(n_samples, n_outputs)

# Create random one-hot labels for classification
for i in 1:n_samples
    Ytrain[i, rand(1:n_outputs)] = 1.0
end

# Dummy test set
Xtest = randn(100, n_features)
Ytest = zeros(100, n_outputs)
for i in 1:100
    Ytest[i, rand(1:n_outputs)] = 1.0
end

network_params = Dict(
    "data" => Dict(
        "Xtrain" => Xtrain,
        "Ytrain" => Ytrain,
        "Xtest"  => Xtest,
        "Ytest"  => Ytest,
        "nFeatures" => n_features,
        "nLabels" => n_outputs
    ),
    "hiddenLayers" => [5, 5],    # example: two hidden layers of 5 neurons
    "HUtype" => "ReLU",
    "OUtype" => "linear"
)

network = Network.Net(network_params)  # Assume dummy placeholder for now

X = randn(n_samples, n_features)
Y = zeros(n_samples, n_outputs)
for i in 1:n_samples
    Y[i, rand(1:n_outputs)] = 1.0  # one-hot encoded
end

neuronsPerLayer = [5, 5, 5, 3]
nLayers = length(neuronsPerLayer)
nParams = sum((neuronsPerLayer[i] + 1) * neuronsPerLayer[i+1] for i in 1:nLayers - 1)
thetavec = randn(nParams)
println(size(thetavec))
designVar = LearnableVars(neuronsPerLayer[1:end-1], nLayers-1, thetavec)  # construct the correct type

# --- Loss Functional ---
lossParams = Dict(
    "costType" => "L2",
    "designVariable" => designVar,
    "network" => network,
    "data" => Dict(
        "Xtrain" => Xtrain,
        "Ytrain" => Ytrain,
        "Xtest" => Xtest,
        "Ytest" => Ytest
    )
)

lossFunc = LossFunctional.LossFunc(lossParams)

# --- L2 Regularization Shape Function ---
shFuncParams = Dict("designVariable" => designVar)
regFunc = Sh_Func_L2norm.ShFuncL2norm(shFuncParams)

# --- CostNN ---
costNNParams = Dict{String, Any}(
    "shapeFunctions" => [lossFunc, regFunc],
    "weights" => [1.0, 0.1]
)

costNN = CostNN.CostNNStruct(costNNParams)

println("Running full batch compute...")
CostNN.computeFunctionAndGradient(costNN, thetavec)
println("Value: ", costNN.value)
#println("Gradient: ", costNN.gradient)

println("\nRunning stochastic compute...")
CostNN.computeStochasticFunctionAndGradient!(costNN, thetavec)
println("Value: ", costNN.value)
println("Gradient: ", costNN.gradient)

# --- Validate early stopping logic ---
alarm = 0.0
minTestError = 1.0
alarm, minTestError = CostNN.validateES(costNN, alarm, minTestError)
println("\nValidateES => Alarm: $alarm | MinTestError: $minTestError")
