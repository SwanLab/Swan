using Random
using LinearAlgebra

# Include your modules (adjust paths if needed)
include("Network/Network.jl")
include("LossFunctional/LossFunctional.jl")
using .Network
using .LossFunctional

# Shortcut for convenience
const Net = Network.Net
const LossFunc = LossFunctional.LossFunc
const computeFunctionAndGradient = LossFunctional.computeFunctionAndGradient
const computeStochasticCostAndGradient = LossFunctional.computeStochasticCostAndGradient

# === Prepare dummy data ===
Random.seed!(1234)

nSamples = 500
nFeatures = 10
nLabels = 3

# Random input matrix: nSamples x nFeatures
Xtrain = randn(nSamples, nFeatures)
Ytrain = zeros(nSamples, nLabels)

# Create random one-hot labels for classification
for i in 1:nSamples
    Ytrain[i, rand(1:nLabels)] = 1.0
end

# Dummy test set
Xtest = randn(100, nFeatures)
Ytest = zeros(100, nLabels)
for i in 1:100
    Ytest[i, rand(1:nLabels)] = 1.0
end

# === Create parameters dictionaries ===

network_params = Dict(
    "data" => Dict(
        "Xtrain" => Xtrain,
        "Ytrain" => Ytrain,
        "Xtest"  => Xtest,
        "Ytest"  => Ytest,
        "nFeatures" => nFeatures,
        "nLabels" => nLabels
    ),
    "hiddenLayers" => [5, 5],    # example: two hidden layers of 5 neurons
    "HUtype" => "ReLU",
    "OUtype" => "softmax"
)

# Create the network object
net = Net(network_params)

# Loss functional params
loss_params = Dict(
    "costType" => "-loglikelihood",  # or "L2"
    "designVariable" => net.learnableVariables,
    "network" => net,
    "data" => network_params["data"]
)

loss = LossFunc(loss_params)

# === Run full batch computation ===
println("\nRunning computeFunctionAndGradient (full batch)...")
cost, grad = computeFunctionAndGradient(loss, net.learnableVariables.thetavec)
println("Cost (full batch): ", cost)
println("Gradient norm (full batch): ", norm(grad))

# === Run stochastic batch computation ===
println("\nRunning computeStochasticCostAndGradient (stochastic batches)...")
for batch_i in 1:5
    cost_s, grad_s, batch_done = computeStochasticCostAndGradient(loss, net.learnableVariables.thetavec, true)
    println("Batch $batch_i cost: ", cost_s)
    println("Batch $batch_i gradient norm: ", norm(grad_s))
    println("Is batch depleted? ", batch_done)
end

# === Check test error ===
test_err = LossFunctional.getTestError(loss)
println("\nTest error: ", test_err)
