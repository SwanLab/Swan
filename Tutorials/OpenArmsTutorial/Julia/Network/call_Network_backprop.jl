using JSON
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .LearnableVariables
using .Network

# Read input
args = JSON.parsefile(ARGS[1])

# Reconstruct Network object from serialized input
net = Network.Net(args["network"])

# Convert inputs
Yb = Matrix{Float64}(args["Yb"])
dLF = Matrix{Float64}(args["dLF"])

# Call method
dc = Network.backprop(net, Yb, dLF)

# Return output
result = Dict("dc" => dc)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
