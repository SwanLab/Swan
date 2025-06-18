using JSON
include("Sh_Func_L2norm.jl")
using .Sh_Func_L2norm

# Read input JSON file passed via command line
args = JSON.parsefile(ARGS[1])

# Construct the object
params = Dict("designVariable" => args["designVariable"])
obj = ShFuncL2norm(params)

# Save the object in a registry (for reuse)
const objectRegistry = Dict{String, ShFuncL2norm}()
objectRegistry[args["id"]] = obj

# Write a confirmation 
open(args["output"], "w") do f
    write(f, JSON.json(Dict("status" => "created", "id" => args["id"])))
end

