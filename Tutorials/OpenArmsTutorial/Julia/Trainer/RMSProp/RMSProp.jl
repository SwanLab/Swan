module RMSProp

export RMSPropStruct, step

using ..SGD  # Import SGD module for SGDStruct

mutable struct RMSPropStruct
    sgd::SGD.SGDStruct
    rho::Float64
    r::Vector{Float64}
end

# Constructor 
function RMSPropStruct(s::Dict{String, Any})
    sgd_obj = SGD.SGDStruct(s)

    # Extract RMSProp-specific hyperparameter
    rho = s["rho"]

    # Initialize r vector to zeros with dimension matching thetavec
    thetavec = s["designVariable"]["thetavec"]
    r = zeros(length(thetavec))

    return RMSPropStruct(sgd_obj, rho, r)
end

# Performs one RMSProp update step
function step(obj::RMSPropStruct, x::Vector{Float64}, e::Float64, grad::Vector{Float64}, F, Xb, Yb)
    # Update running average of squared gradients
    obj.r = obj.rho .* obj.r .+ (1 .- obj.rho) .* (grad .^ 2)

    # Compute parameter update with RMSProp scaling
    x = x .- e ./ sqrt.(obj.r .+ 1e-6) .* grad

    return x, grad
end

end