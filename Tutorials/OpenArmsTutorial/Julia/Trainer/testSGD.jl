# --- TestSGD.jl ---

using Random

# Load your modules (adjust paths if needed)
include("Trainer.jl")
include("SGD/SGD.jl")

using .Trainer
using .SGD
using LinearAlgebra

# Mock Objective Function type with minimal required API
mutable struct MockObjectiveFunction
    thetavec::Vector{Float64}
    # We add dummy values for attributes that might be accessed in Trainer/SGD
    regularization::Float64
    loss::Float64
    value::Float64
    gradient::Vector{Float64}
    batchDepleted::Bool
end

function MockObjectiveFunction()
    # start with zero vector
    theta = randn(5)
    return MockObjectiveFunction(theta, 0.0, 0.0, 0.0, zeros(length(theta)), false, )
end

# Define minimal interface used by SGD (simulate stochastic function & gradient)
function computeStochasticFunctionAndGradient(obj::MockObjectiveFunction, x::Vector{Float64})
    # Example: simple quadratic function and gradient
    obj.value = sum(x.^2)
    obj.gradient = 2 .* x
    # Simulate batch depletion toggle for demo
    obj.batchDepleted = !obj.batchDepleted
    return obj.value, obj.gradient
end

# Set batch depletion status
function isBatchDepleted(obj::MockObjectiveFunction)
    return obj.batchDepleted
end

# For validation early stopping dummy function
function validateES(obj::MockObjectiveFunction, alarm::Int, minTestError::Float64)
    return alarm, minTestError
end

# Create a dictionary mimicking your parameters for SGD
params = Dict(
    "costFunc" => Dict(
        "thetavec" => randn(5),
        "computeStochasticFunctionAndGradient" => (x->computeStochasticFunctionAndGradient(x[1], x[2])), # We'll adapt this below
        "isBatchDepleted" => false,
        "validateES" => validateES,
        "setBatchMover" => (moveBatch->nothing),
        "value" => 0.0,
        "gradient" => zeros(5),
        "regularization" => 0.0,
        "loss" => 0.0
    ),
    "designVariable" => Dict(
        "thetavec" => randn(5)
    ),
    "maxEpochs" => 3,
    "learningRate" => 0.1,
    "plotter" => nothing
)

# Fix the function in costFunc dict for compatibility with SGD call:
# We'll wrap the MockObjectiveFunction itself for the call

mockObj = MockObjectiveFunction()
params["costFunc"] = mockObj

# Create the SGD struct from params
sgd = SGD.SGDStruct(params)

# Patch the computeStochasticFunctionAndGradient to delegate correctly:
function SGD.computeStochasticFunctionAndGradient(obj::SGD.SGDStruct, theta::Vector{Float64}, moveBatch::Bool)
    # update batch depletion for test
    obj.trainer.objectiveFunction.batchDepleted = isBatchDepleted(obj.trainer.objectiveFunction)
    # call the mock compute function
    f, grad = computeStochasticFunctionAndGradient(obj.trainer.objectiveFunction, theta)
    return f, grad
end

# Also patch isBatchDepleted check inside SGD.optimize loop:
function SGD.optimize(obj::SGD.SGDStruct, th0::Vector{Float64})
    epsilon      = obj.learningRate
    iter         = -1
    funcount     = 0
    alarm        = 0
    minTestError = 1.0

    KPI = Dict(
        :epoch => 1,
        :alarm => 0,
        :gnorm => 1.0,
        :cost  => 1.0
    )

    theta = th0
    while !SGD.isCriteriaMet(obj, KPI)
        state = iter == -1 ? :init : :iter
        newEpoch = true
        moveBatch = true

        while !obj.trainer.objectiveFunction.batchDepleted || newEpoch
            f, grad = SGD.computeStochasticFunctionAndGradient(obj, theta, moveBatch)
            epsilon, theta, funcount = SGD.lineSearch(obj, theta, grad, f, epsilon, funcount)

            funcount += 1
            iter += 1
            newEpoch = false

            KPI[:cost]  = f
            KPI[:gnorm] = norm(grad)
            SGD.displayIter(obj, iter, funcount, theta, epsilon, state, KPI)
        end

        KPI[:epoch] += 1
        kpi_alarm, minTestError = validateES(obj.trainer.objectiveFunction, alarm, minTestError)
        KPI[:alarm] = kpi_alarm
    end
end

println("Running SGD test optimization...")
SGD.compute(sgd)
println("SGD test run complete.")
