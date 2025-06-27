using JSON
include("Trainer.jl")
include("SGD/SGD.jl")        # Include whichever trainers you support
include("Fminunc/Fminunc.jl")
include("Nesterov/Nesterov.jl")
include("RMSProp/RMSProp.jl")

using .Trainer
using .SGD       # Assuming your SGD is in module SGDModule
using .Fminunc
using .Nesterov
using .RMSProp

# Step 1: Read JSON input
args = JSON.parsefile(ARGS[1])

# Step 2: Create Trainer object using factory function
trainer = Trainer.Create(args)

# Step 3: Prepare output
result = Dict(
    "status" => "created",
    "type" => args["type"]
)

# Step 4: Write to JSON output file
open(args["output"], "w") do f
    write(f, JSON.json(result))
end
