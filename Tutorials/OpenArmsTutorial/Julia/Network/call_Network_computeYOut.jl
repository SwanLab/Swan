using JSON
include("../LearnableVariables/LearnableVariables.jl")  # Adjust if needed
include("Network.jl")
using .Network
using .LearnableVariables

# Step 1: Read JSON input
args = JSON.parsefile(ARGS[1])

# Step 2: Reconstruct the Net object
net = Network.Net(args)

# Step 2.5: Extract input data
Xb = hcat([Float64.(row) for row in args["Xb"]]...) 

# Step 3: Call method
yOut = Network.computeYOut(net, Xb)

println("I called computeYOut")

# Step 4: Write result into output file
result = Dict("yOut" => yOut)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
