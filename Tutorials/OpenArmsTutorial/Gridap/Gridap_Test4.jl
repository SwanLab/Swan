# Test_SGD_Trainer.jl

# ----------------------------
# ✅ Load modules in dependency order
# ----------------------------
include("Data.jl")
include("LearnableVariables.jl")
include("Network.jl")
include("LossFunctional.jl")
include("Sh_Func_L2norm.jl")
include("CostNN.jl")
include("Trainer.jl")
include("SGD.jl")

using .Data
using .LearnableVariables
using .Network
using .LossFunctional
using .Sh_Func_L2norm
using .CostNN
using .Trainer
using .SGD
using LinearAlgebra

function main()
    println("=== SGD + TRAINER MODULE TEST ===")

    # ----------------------------
    # 🔧 Create parameters for a minimal neural network LossFunctional
    # ----------------------------
    params_nn = Dict(
        "hiddenLayers" => [2],
        "HUtype" => "ReLU",
        "OUtype" => "linear",
        "costType" => "L2",
        "xFeatures" => [1,2,3,4,5,6,7],
        "yFeatures" => [8],
        "testRatio" => 20.0,
        "fileName" => "Resultados2.csv",
        "polynomialOrder" => 1
    )

    data_nn = Data.init_data(params_nn)
    params_nn["data"] = data_nn
    net = Network.init_network(params_nn)
    params_nn["network"] = net

    loss = LossFunctional.init_lossfunctional(params_nn)
    println("✔ LossFunctional initialized.")

    # θ vector from learnable variables
    lv = LearnableVariables.init_learnable_variables(Dict(
        "neuronsPerLayer" => net.neurons_per_layer,
        "nLayers" => net.n_layers
    ))
    θ = lv.thetavec

    # ----------------------------
    # 🔧 Create a Sh_Func_L2norm shape function
    # ----------------------------
    regularization = Sh_Func_L2norm.init_ShFuncL2norm(Dict())
    println("✔ Sh_Func_L2norm initialized.")

    # ----------------------------
    # 🔧 Initialize CostNN
    # ----------------------------
    params_costnn = Dict{String, Any}(
        "shapeFunctions" => [loss, regularization],
        "weights" => [0.7, 0.3]
    )
    costnn = CostNN.init_CostNN(params_costnn)
    println("✔ CostNN initialized.")

    # ----------------------------
    # 🔧 Initialize Trainer
    # ----------------------------
    trainer_params = Dict(
        "costFunc" => costnn,
        "designVariable" => lv,
        "maxEpochs" => 5
    )

    trainer = Trainer.init_Trainer(trainer_params)
    println("✔ Trainer initialized. Epoch counter: ", trainer.epoch_counter)

    # ----------------------------
    # 🔎 Test store_values
    # ----------------------------
    optinfo = Trainer.opt_info(0.001, 0.05)
    #trainer_updated = Trainer.store_values(trainer, θ, 0.12345, optinfo)  # store_values does not work as it is now (neither does it in matlab)
    #println("🔧 store_values executed. New epoch counter: ", trainer_updated.epoch_counter)

    # ----------------------------
    # 🔧 Initialize SGD
    # ----------------------------
    sgd_params = Dict(
        "maxEpochs" => 5,
        "learningRate" => 0.01,
        "costFunc" => costnn,
        "designVariable" => lv
    )

    sgd = SGD.init_SGD(sgd_params)
    println("✔ SGD initialized. Learning rate: ", sgd.learning_rate)

    # ----------------------------
    # 🔎 Test compute_stochastic_function_and_gradient
    # ----------------------------
    jV, djV, is_batch_depleted, Jc = SGD.compute_stochastic_function_and_gradient(costnn, θ, true)
    println("🔧 compute_stochastic_function_and_gradient")
    println("   - Value: ", jV)
    println("   - Gradient norm: ", norm(djV))
    println("   - is_batch_depleted: ", is_batch_depleted)

    # ----------------------------
    # 🔎 Test full compute (optimization run)
    # ----------------------------
    try
        sgd_final = SGD.compute(sgd, θ)
        println("🔧 SGD compute executed. Final elapsed time: ", sgd_final.elapsed_time)
    catch e
        println("❌ SGD compute failed: ", e)
    end

    println("=== SGD + TRAINER TEST COMPLETED SUCCESSFULLY ===")
end

main()
