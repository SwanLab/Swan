using JSON
include("Data.jl")
using .Data

function matrix_to_nested_array(mat)
    s = size(mat,1)
    return [vec(mat[i, :]) for i in 1:s]
end

function vector_to_array(vec)
    return collect(vec)
end

# Step 1: Read input JSON file
args = JSON.parsefile(ARGS[1])

if isa(args["yFeatures"], Number)
    args["yFeatures"] = [args["yFeatures"]]  # convert scalar to 1-element vector
end

# Step 2: Construct Data object
params = args
obj = Data.DataStruct(params)
println(size(obj.Xtrain))
# Step 3: Prepare output dictionary with all public properties converted to JSON arrays
result = Dict(
    "status" => "created",
    "nFeatures" => obj.nFeatures,
    "nSamples" => obj.nSamples,
    "nLabels" => obj.nLabels,
    "Xtrain" => matrix_to_nested_array(obj.Xtrain),
    "Ytrain" => matrix_to_nested_array(obj.Ytrain),
    "Xtest"  => matrix_to_nested_array(obj.Xtest),
    "Ytest"  => matrix_to_nested_array(obj.Ytest),
    "Ntest"  => obj.Ntest,
    "batchSize" => obj.batchSize,
    "Batch_nD"  => obj.Batch_nD,
    "Batch_nB"  => obj.Batch_nB,
    "muX"      => vector_to_array(obj.muX),
    "sigmaX"   => vector_to_array(obj.sigmaX),
    "muY"      => vector_to_array(obj.muY),
    "sigmaY"   => vector_to_array(obj.sigmaY)
)

# Step 4: Write JSON output
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
