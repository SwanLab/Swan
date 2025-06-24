using JSON
include("../LearnableVariables/LearnableVariables.jl")
include("Network.jl")
using .LearnableVariables
using .Network

# Read input
args = JSON.parsefile(ARGS[1])

# Reconstruct Network object from serialized input
net = Network.Net(args)

# Convert inputs
Xb = hcat(args["Xb"]...)                # <--- Required input to computeAvalues
Xb = reshape(Xb, :, net.nPolyFeatures) 

Yb = hcat(args["Yb"]...)
Yb = reshape(Yb, :, 1)
println("typeof(args[\"dLF\"]): ", typeof(args["dLF"]))
println("size(args[\"dLF\"]): ", size(args["dLF"]))
dLF = hcat(args["dLF"]...)
dLF = reshape(dLF, :, 1)

# Compute aValues first 
Network.computeAvalues(net, Xb)

# Call method
dc = Network.backprop(net, Yb, dLF)

# Return output
result = Dict("dc" => dc)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
