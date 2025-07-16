module MyProject

include("LearnableVariables.jl")
using .LearnableVariables
export LearnableVariables

include("Data.jl")
using .Data
export Data

include("Network.jl")
using .Network
export Network

include("LossFunctional.jl")
using .LossFunctional
export LossFunctional

include("Sh_Func_L2norm.jl")
using .Sh_Func_L2norm
export Sh_Func_L2norm

include("CostNN.jl")
using .CostNN
export CostNN

include("PlotterNN.jl")
using .PlotterNN
export PlotterNN

include("Trainer.jl")
using .Trainer
export Trainer

include("SGD.jl")
using .SGD
export SGD

include("Nesterov.jl")
using .Nesterov
export Nesterov

include("RMSProp.jl")
using .RMSProp
export RMSProp

include("Fminunc.jl")
using .Fminunc
export Fminunc

include("OptimizationProblemNN.jl")
using .OptimizationProblemNN
export OptimizationProblemNN

end
