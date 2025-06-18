using MATLAB

# Initialization of hyperparameters
pol_deg         = 1
testratio       = 20
lambda          = 0.0
learningRate    = 0.2
hiddenLayers    = fill(128, 6)

# INITIALIZATION 
s = Dict{Symbol, Any}()

# Store dataset file name
s[:fileName] = "Tutorials/OpenArmsTutorial/Julia/Resultados2.csv"

# Load model parameters
s[:polynomialOrder] = pol_deg
s[:testRatio]        = testratio

s[:networkParams] = Dict(
    :hiddenLayers => hiddenLayers,
    :HUtype => "ReLU",
    :OUtype => "linear"
)

s[:optimizerParams] = Dict(
    :learningRate => learningRate,
    :maxEpochs    => 3000
)

s[:costParams] = Dict(
    :lambda   => lambda,
    :costType => "L2"
)
# Select the model's features
s[:xFeatures] = [1, 2, 3, 4, 5, 6, 7]
s[:yFeatures] = [8]



mat"""
addpath('C:/Users/jsend/OneDrive/Documentos/GitHub/Swan');
"""
# ========== Call Data class in MATLAB ==========
mat"data = Data(s);"
mat"s.data = data;"

# ========== Train the model ==========
mat"opt = OptimizationProblemNN(s);"
mat"opt.solve();"
mat"opt.plotCostFnc();"

# ========== Get test data ==========
mat"[Xtest, Ytest] = opt.getTestData();"

# ========== Predict using the network ==========
mat"""
    network = opt.getNetwork();
    Ypred = zeros(size(Ytest));
    for i = 1:size(Xtest, 1)
        Ypred(i) = network.computeYOut(Xtest(i, :));
    end
"""

# ========== Histogram Plots ==========
Ypred = mat"Ypred"
Ytest = mat"Ytest"
histogram(Ypred, bins=30, xlims=(-1, 2), title="Distribution of predicted Ytest")
histogram(Ytest, bins=30, xlims=(-1, 2), title="Distribution of Test Y")

# ========== Denormalization ==========
mat"""
    Xtest = Xtest .* data.sigmaX + data.muX;
    Ypred = Ypred .* data.sigmaY + data.muY;
    Ytest = Ytest .* data.sigmaY + data.muY;
    difference = Ytest - Ypred;
"""

# ========== Plot: Fuel consumption vs speed^3 ==========
Xtest = mat"Xtest"
Ytest = mat"Ytest"
plot(Xtest[:, 4], Ytest, seriestype=:scatter,
     xlabel="Speed cubed (m/s)^3", ylabel="Fuel consumption")

# ========== Mean Squared Error ==========
Ypred = mat"Ypred"
mse = mean((Ypred .- Ytest).^2)
println("Error cuadrático medio (MSE) en los datos de prueba: $mse")

# ========== Build output DataFrame ==========
difference = mat"difference"

# Use original Xtest values for input table
input_data = DataFrame(Xtest, :auto)
rename!(input_data, [:rpm, "Windy(cosine)", "Windy(m/s)", "Speed^3 (m/s)^3", :Yaw, :Pitch, :Roll])

output_data = DataFrame(
    "Cons. real"       => Ytest,
    "Cons. prediction" => Ypred,
    "Difference"       => difference
)

result_table = hcat(input_data, output_data)

println("\nTabla de resultados:")
show(result_table, allcols=true)


