using JSON
include("Trainer.jl")
using .Trainer

# Step 1: Read JSON input
args = JSON.parsefile(ARGS[1])

# Step 2: Construct TrainerStruct
objFunc = args["objectiveFunction"]
designVar = args["designVariable"]
trainer = Trainer.TrainerStruct(objFunc, designVar)

# Step 3: Return confirmation and empty histories for debugging or inspection
result = Dict(
    "status" => "created",
    "costHist" => trainer.costHist,  # should be empty (0,3) matrix
    "optHist" => trainer.optHist     # should be empty (0,2) matrix
)

# Step 4: Write output
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
