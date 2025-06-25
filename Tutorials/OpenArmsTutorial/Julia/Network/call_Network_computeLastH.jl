using JSON
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .Network
using .LearnableVariables

# Parse input arguments
args = JSON.parsefile(ARGS[1])

# Reconstruct Net object from JSON
net = Net(args)

# Load input matrix X
X = Matrix(hcat(map(x -> Float64.(x), args["X"])...)')

# Compute the last hidden layer activation
H = computeLastH(net, X)

println("I called computeLastH")
# Write result
result = Dict("H" => H)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
