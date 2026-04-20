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

include("Adam.jl")
using .Adam
export Adam

include("OptimizationProblemNN.jl")
using .OptimizationProblemNN
export OptimizationProblemNN

end