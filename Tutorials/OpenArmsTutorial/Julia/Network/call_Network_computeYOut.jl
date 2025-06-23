using JSON
include("../LearnableVariables/LearnableVariables.jl")  # Adjust if needed
include("Network.jl")
using .Network
using .LearnableVariables

# Read JSON input
args = JSON.parsefile(ARGS[1])

# Reconstruct the Net object
net = Net(args)

# Extract input data
Xb = hcat([Float64.(row) for row in args["Xb"]]...)  # Transpose if needed

# Call method
yOut = computeYOut(net, Xb)

# Write result
result = Dict("yOut" => yOut)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
