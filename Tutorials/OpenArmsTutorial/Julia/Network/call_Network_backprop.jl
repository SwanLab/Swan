using JSON
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .LearnableVariables
using .Network

# Step 1: Read JSON input
args = JSON.parsefile(ARGS[1])

# Step 2: Reconstruct the Net object
net = Network.Net(args)

# Step 2.5: Extract input data
Xb = hcat(args["Xb"]...)                # <--- Required input to computeAvalues
Xb = reshape(Xb, :, net.nPolyFeatures) 

Yb = hcat(args["Yb"]...)
Yb = reshape(Yb, :, 1)
dLF = hcat(args["dLF"])

# Step 3: Call method

# Compute aValues first 
Network.computeAvalues(net, Xb)

dc = Network.backprop(net, Yb, dLF)

println("I called backprop")

# Step 4: Write result into output file
result = Dict("dc" => dc)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
