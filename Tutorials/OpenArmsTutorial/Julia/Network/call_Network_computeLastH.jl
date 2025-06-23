using JSON
include("Network.jl")
include("../LearnableVariables/LearnableVariables.jl")
using .Network
using .LearnableVariables

# Parse input arguments
args = JSON.parsefile(ARGS[1])

# Reconstruct Net object from JSON
net = Net(args)

# Load input matrix X
X = Matrix{Float64}(args["X"])

# Compute the last hidden layer activation
H = computeLastH(net, X)

# Write result
result = Dict("H" => H)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
