using JSON
include("Sh_Func_L2norm.jl")
using .Sh_Func_L2norm

# Step 1: Read the JSON input file
args = JSON.parsefile(ARGS[1])

# Step 2: Construct the object
thetavec = Vector{Float64}(args["designVariable"]["thetavec"])
designVariable = Dict("thetavec" => thetavec)
params = designVariable
obj = Sh_Func_L2norm.ShFuncL2norm(params)

# Step 3: Return basic confirmation
result = Dict(
    "status" => "created",
    "thetavec" => obj.designVariable["thetavec"]
)

# Step 4: Write to output JSON
open(args["output"], "w") do f
    write(f, JSON.json(result))
end

