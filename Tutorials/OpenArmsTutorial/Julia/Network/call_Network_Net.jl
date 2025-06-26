using JSON
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .LearnableVariables
using .Network


# Read input
args = JSON.parsefile(ARGS[1])

# Create Julia object
net = Network.Net(args)

# Prepare the minimal data to return to MATLAB
# (MATLAB doesn't need to hold the full object, they will be defined through other methods)
result = Dict(
    "neuronsPerLayer" => net.neuronsPerLayer
)

# Write output
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
