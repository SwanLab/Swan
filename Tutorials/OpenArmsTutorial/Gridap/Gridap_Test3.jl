# Test_CostNN.jl

# ----------------------------
# ✅ Load modules in dependency order
# ----------------------------
include("Data.jl")
include("LearnableVariables.jl")
include("Network.jl")
include("LossFunctional.jl")
include("Sh_Func_L2norm.jl")
include("CostNN.jl")

using .Data
using .LearnableVariables
using .Network
using .LossFunctional
using .Sh_Func_L2norm
using .CostNN
using LinearAlgebra

function main()
    println("=== COSTNN MODULE TEST ===")

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

    lf = LossFunctional.init_lossfunctional(params_nn)

    # θ vector from learnable variables
    lv = LearnableVariables.init_learnable_variables(Dict(
        "neuronsPerLayer" => net.neurons_per_layer,
        "nLayers" => net.n_layers
    ))
    θ = lv.thetavec

    println("✔ LossFunctional initialized.")

    # ----------------------------
    # 🔧 Create a Sh_Func_L2norm shape function
    # ----------------------------
    #shparams = Dict("designVariable" => Dict("thetavec" => θ))
    shfunc = Sh_Func_L2norm.init_ShFuncL2norm(Dict())
    println("✔ Sh_Func_L2norm initialized.")

    # ----------------------------
    # 🔧 Initialize CostNN with both shape functions
    # ----------------------------
    params_costnn = Dict{String, Any}(
        "shapeFunctions" => [lf, shfunc],
        "weights" => [0.7, 0.3]
    )
    costnn = CostNN.init_CostNN(params_costnn)
    println("✔ CostNN initialized.")

    # ----------------------------
    # 🔎 Test compute_function_and_gradient
    # ----------------------------
    jV, djV, Jc = CostNN.compute_function_and_gradient(costnn, θ)
    println("🔧 compute_function_and_gradient")
    println("   - Value: ", jV)
    println("   - Gradient norm: ", norm(djV))

    # ----------------------------
    # 🔎 Test compute_stochastic_function_and_gradient
    # ----------------------------
    jVs, djVs, isBD, Jcs = CostNN.compute_stochastic_function_and_gradient(costnn, θ)
    println("🔧 compute_stochastic_function_and_gradient")
    println("   - Value: ", jVs)
    println("   - Gradient norm: ", norm(djVs))
    println("   - isBatchDepleted: ", isBD)

    # ----------------------------
    # 🔎 Test obtain_number_fields
    # ----------------------------
    n_fields = CostNN.obtain_number_fields(costnn)
    println("🔧 obtain_number_fields: ", n_fields)

    # ----------------------------
    # 🔎 Test get_title_fields
    # ----------------------------
    titles = CostNN.get_title_fields(costnn)
    println("🔧 get_title_fields: ", titles)

    # ----------------------------
    # 🔎 Test get_fields
    # ----------------------------
    field1 = CostNN.get_fields(Jc, 1)
    println("🔧 get_fields[1]: ", field1)

    # ----------------------------
    # 🔎 Test set_batch_mover
    # ----------------------------
    costnn2 = CostNN.set_batch_mover(costnn, false)
    println("🔧 set_batch_mover: moveBatch changed to ", costnn2.moveBatch)

    # ----------------------------
    # 🔎 Test validate_ES
    # ----------------------------
    alarm = 0.0
    minTestError = 1.0
    alarm_new, minTestError_new = CostNN.validate_ES(costnn, alarm, minTestError, θ)
    println("🔧 validate_ES")
    println("   - New alarm: ", alarm_new)
    println("   - New minTestError: ", minTestError_new)

    println("=== COSTNN TEST COMPLETED SUCCESSFULLY ===")
end

main()
