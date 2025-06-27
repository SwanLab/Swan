using JSON
include("SGD.jl")
using .SGD

# Step 1: Read JSON input (here, if needed)
args = JSON.parsefile(ARGS[1])

# Step 2: Reconstruct SGD object if necessary, or use a global if your setup maintains it in the same session
obj = SGD.SGDStruct(args)

# Step 3: Call compute
SGD.compute(obj)

# Step 4: Return simple status
result = Dict("status" => "compute done")
open(args["output"], "w") do f
    write(f, JSON.json(result))
end