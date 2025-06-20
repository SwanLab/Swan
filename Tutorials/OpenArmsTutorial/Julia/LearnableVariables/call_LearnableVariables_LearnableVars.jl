using JSON
include("LearnableVariables.jl")
using .LearnableVariables

# Step 1: Read JSON input
args = JSON.parsefile(ARGS[1])

# Step 2: Construct LearnableVars object
params = args
obj = LearnableVariables.LearnableVars(params)

# Step 3: Write output
result = Dict(
    "status" => "created",
    "thetavec" => obj.thetavec
)

open(args["output"], "w") do f
    write(f, JSON.json(result))
end
