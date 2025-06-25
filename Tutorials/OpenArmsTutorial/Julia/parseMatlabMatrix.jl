function parseMatlabMatrix(arr::Vector{Any})
    return hcat(map(x -> Float64.(x), arr)...)'  # Convert and transpose into M×N
end