using JSON
include("Network.jl")
include("../LearnableVariables/LearnableVariables.jl")
using .Network
using .LearnableVariables

# Read input from JSON
args = JSON.parsefile(ARGS[1])

# Reconstruct Net object
net = Net(args)

# Get LearnableVars object
lvars = getLearnableVariables(net)

# Prepare result as a Dict for serialization
result = Dict(
    "neuronsPerLayer" => lvars.neuronsPerLayer,
    "nLayers" => lvars.nLayers,
    "thetavec" => lvars.thetavec
)

# Write to output file
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
