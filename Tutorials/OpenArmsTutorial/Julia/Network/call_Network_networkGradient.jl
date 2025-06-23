using JSON
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .Network
using .LearnableVariables

# Parse input arguments
args = JSON.parsefile(ARGS[1])

# Reconstruct Net object from input parameters
net = Net(args)

# Convert input data
X = Matrix{Float64}(args["X"])

# Call the networkGradient method
dy = networkGradient(net, X)

# Write output to JSON file
result = Dict("dy" => dy)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
