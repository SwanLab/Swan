using JSON
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .Network
using .LearnableVariables

# Step 1: Read JSON input
args = JSON.parsefile(ARGS[1])

# Step 2: Reconstruct the Net object
net = Network.Net(args)

# Step 3: Call method
lvars = Network.getLearnableVariables(net)

# Step 3.5: Prepare result as a Dict for serialization
result = Dict(
    "neuronsPerLayer" => lvars.neuronsPerLayer,
    "nLayers" => lvars.nLayers,
    "thetavec" => lvars.thetavec
)
println("I called getLearnableVariables")

# Step 4: Write result into output file
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
