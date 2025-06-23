using JSON
include("Network.jl")
include("../LearnableVariables/LearnableVariables.jl")  # Adjust if needed
using .Network
using .LearnableVariables

# Read JSON input
args = JSON.parsefile(ARGS[1])

# Reconstruct the Net object
net = Net(args["network"])

# Extract input data
Xb = Matrix{Float64}(args["Xb"])

# Call method
yOut = computeYOut(net, Xb)

# Write result
result = Dict("yOut" => yOut)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
