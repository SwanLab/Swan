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
network = nothing  # Assume dummy placeholder for now

n_features = 5
n_outputs = 3
n_samples = 20

X = randn(n_samples, n_features)
Y = zeros(n_samples, n_outputs)
for i in 1:n_samples
    Y[i, rand(1:n_outputs)] = 1.0  # one-hot encoded
end

neuronsPerLayer = [5, 5, 5]
nLayers = length(neuronsPerLayer)
thetavec = randn(n_features * n_outputs)
designVar = LearnableVars(neuronsPerLayer, nLayers, thetavec)  # construct the correct type

# --- Loss Functional ---
lossParams = Dict(
    "costType" => "L2",
    "designVariable" => designVar,
    "network" => network,
    "data" => Dict(
        "Xtrain" => X,
        "Ytrain" => Y,
        "Xtest" => X,
        "Ytest" => Y
    )
)

lossFunc = LossFunctional.LossFunc(lossParams)

# --- L2 Regularization Shape Function ---
shFuncParams = Dict("designVariable" => designVar)
regFunc = Sh_Func_L2norm.ShFuncL2norm(shFuncParams)

# --- CostNN ---
costNNParams = Dict(
    "shapeFunctions" => [lossFunc, regFunc],
    "weights" => [1.0, 0.1]
)

costNN = CostNN.CostNNStruct(costNNParams)

# Initial x
x = randn(n_features * n_outputs)

println("Running full batch compute...")
CostNN.computeFunctionAndGradient(costNN, x)
println("Value: ", costNN.value)
println("Gradient: ", costNN.gradient)

println("\nRunning stochastic compute...")
CostNN.computeStochasticFunctionAndGradient!(costNN, x)
println("Value: ", costNN.value)
println("Gradient: ", costNN.gradient)

# --- Validate early stopping logic ---
alarm = 0.0
minTestError = 1.0
alarm, minTestError = CostNN.validateES(costNN, alarm, minTestError)
println("\nValidateES => Alarm: $alarm | MinTestError: $minTestError")
