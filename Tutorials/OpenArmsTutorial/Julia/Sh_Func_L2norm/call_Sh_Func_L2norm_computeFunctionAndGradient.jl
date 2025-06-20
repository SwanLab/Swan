using JSON
include("Sh_Func_L2norm.jl")
using .Sh_Func_L2norm

# Read input
args = JSON.parsefile(ARGS[1])

# Create object
obj = Sh_Func_L2norm.ShFuncL2norm(Dict("designVariable" => args["designVariable"]))

# Convert to Float64
x = Vector{Float64}(args["x"])

# Call method
j, dj = computeFunctionAndGradient(obj, x)

# Write result
result = Dict("j" => j, "dj" => dj)
open(args["output"], "w") do f
    write(f, JSON.json(result))
end


