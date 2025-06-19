using JSON
include("Sh_Func_L2norm.jl")
using .Sh_Func_L2norm

# Read input
args = JSON.parsefile(ARGS[1])

# Create object
obj = Sh_Func_L2norm.ShFuncL2norm(Dict("designVariable" => args["designVariable"]))

# Call method
j, dj = computeFunctionAndGradient(obj, args["x"])

# Write result
result = Dict("j" => j, "dj" => dj)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end