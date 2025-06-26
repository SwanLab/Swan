using JSON

# Include and import modules
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .Network
using .LearnableVariables

# Read input from JSON
args = JSON.parsefile(ARGS[1])

# Reconstruct Net object from JSON arguments
net = Net(args)

# Extract thetavec from args
thetavec = args["thetavec"]

# Call the target method
setLearnableVariables!(net, thetavec)

# Optional: print confirmation
println("Called setLearnableVariables! on net")

# Return empty result (can be expanded if needed)
result = Dict()

# Write result to output JSON file
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
