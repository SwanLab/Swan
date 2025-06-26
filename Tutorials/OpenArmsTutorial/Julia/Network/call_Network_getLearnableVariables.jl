using JSON
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .Network
using .LearnableVariables

# Read input from JSON
args = JSON.parsefile(ARGS[1])

# Reconstruct Net object
net = Network.Net(args)

# Get LearnableVars object
lvars = Network.getLearnableVariables(net)

# Prepare result as a Dict for serialization
result = Dict(
    "neuronsPerLayer" => lvars.neuronsPerLayer,
    "nLayers" => lvars.nLayers,
    "thetavec" => lvars.thetavec
)
println("I called getLearnableVariables")
# Write to output file
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
