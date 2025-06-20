using JSON
include("LearnableVariables.jl")
using .LearnableVariables

# Step 1: Read JSON input
args = JSON.parsefile(ARGS[1])

# Step 2: Reconstruct LearnableVars object
thetavec = Vector{Float64}(args["thetavec"])
nPL = Vector{Int}(args["neuronsPerLayer"])
nLayers = args["nLayers"]
obj = LearnableVariables.LearnableVars(nPL, nLayers, thetavec)

# Step 3: Call the method
W, b = LearnableVariables.reshapeInLayerForm(obj)

# Step 4: Convert output to JSON-serializable format
# Matrices must be converted to nested arrays explicitly
W_json = [collect(eachcol(w)) for w in W]  # column-wise array of vectors
b_json = b  # already a vector of vectors

result = Dict(
    "W" => W_json,
    "b" => b_json
)

open(args["output"], "w") do f
    write(f, JSON.json(result))
end
