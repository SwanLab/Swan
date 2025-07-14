# Gridap_Test2.jl

# ----------------------------
# ✅ Load modules in dependency order
# ----------------------------
include("Data.jl")
include("LearnableVariables.jl")
include("Network.jl")
include("LossFunctional.jl")

using .Data
using .LearnableVariables
using .Network
using .LossFunctional
using LinearAlgebra
using Random

# ----------------------------
# ✅ Define main testing function
# ----------------------------
function main()
    println("=== GRIDAP MIGRATION TEST 1 ===")

    # ------------------------
    # 🔧 PARAMETERS
    # ------------------------
    params = Dict(
        "fileName"         => "Resultados2.csv",
        "testRatio"        => 30.0,
        "polynomialOrder"  => 1,
        "xFeatures"        => [1, 2, 3, 4, 5, 6, 7],
        "yFeatures"        => [8],
        "HUtype"           => "ReLU",
        "OUtype"           => "linear",
        #"neuronsPerLayer"  => [5, 10, 1],   # Example: 5 inputs → 10 hidden → 1 output
        #"nLayers"          => 3,
        "hiddenLayers"     => [10],
        "costType"         => "L2"
    )

    # ------------------------
    # 🧩 INITIALIZE DATA
    # ------------------------
    data = Data.init_data(params)
    println("✔ Data initialized.")
    println("   - Xtrain size: ", size(data.Xtrain))
    println("   - Ytrain size: ", size(data.Ytrain))
    params["data"] = data


    # ------------------------
    # 🧩 INITIALIZE NETWORK
    # ------------------------
    net = Network.init_network(params)
    println("✔ Network initialized.")
    println("   - Layers: ", net.neurons_per_layer)
    params["network"] = net

    # ------------------------
    # 🧩 INITIALIZE LEARNABLE VARIABLES
    # ------------------------
    lv = Network.get_learnable_variables(params["network"])
    println("✔ Learnable variables initialized.")
    println("   - θ vector length: ", length(lv.thetavec))
    params["learnableVariables"] = lv

    
    # ------------------------
    # 🧩 INITIALIZE LOSS FUNCTIONAL
    # ------------------------
    lf = LossFunctional.init_lossfunctional(params)
    println("✔ Loss functional initialized with cost type: ", lf.cost_type)

    # ------------------------
    # 🔎 TEST FORWARD + BACKWARD
    # ------------------------
    θ = lv.thetavec
    j, grad = LossFunctional.compute_loss_and_gradient(lf, θ)
    println("🔧 Full-batch loss computed:")
    println("   - Loss value: ", j)
    println("   - Gradient norm: ", norm(grad))

    # ------------------------
    # 🔎 TEST MINI-BATCH COMPUTATION
    # ------------------------
    order = randperm(size(data.Xtrain, 1))
    i_batch = 1
    j_stoch, grad_stoch, is_depleted = LossFunctional.compute_stochastic_loss_and_gradient(lf, θ, i_batch, order, true)
    println("🔧 Mini-batch loss computed:")
    println("   - Loss value: ", j_stoch)
    println("   - Gradient norm: ", norm(grad_stoch))
    println("   - Dataset depleted: ", is_depleted)

    # ------------------------
    # 🔎 TEST TEST-ERROR (classification tasks)
    # ------------------------
    err = LossFunctional.get_test_error(lf, θ)
    println("✔ Test error: ", err)

    println("=== TEST COMPLETED SUCCESSFULLY ===")
end

# ----------------------------
# ✅ Run main when script executed
# ----------------------------
main()
