module Data

export DataStruct, init_data, plot_data, plot_corr_row, plot_corr_matrix, update_hyperparameter!

using CSV, DelimitedFiles, Statistics, Random, LinearAlgebra

struct DataStruct
    n_features::Int
    n_samples::Int
    n_labels::Int
    Xtrain::Matrix{Float64}
    Ytrain::Matrix{Float64}
    Xtest::Matrix{Float64}
    Ytest::Matrix{Float64}
    Ntest::Int
    batch_size::Int
    batch_nD::Int
    batch_nB::Int
    muX::Matrix{Float64}
    sigmaX::Matrix{Float64}
    muY::Matrix{Float64}
    sigmaY::Matrix{Float64}
    X::Matrix{Float64}
    Y::Matrix{Float64}
    polynomial_order::Int
    data::Matrix{Float64}
    file_name::String
    test_ratio::Float64
    x_features::Vector{Int}
    y_features::Vector{Int}
end

function init_data(cparams::Dict{String, Any})
    file_name        = cparams["fileName"]
    test_ratio       = cparams["testRatio"]
    polynomial_order = cparams["polynomialOrder"]
    x_features       = cparams["xFeatures"]
    y_features       = cparams["yFeatures"]

    data = readdlm(file_name, ',', Float64)
    x = data[:, x_features]
    y = data[:, y_features]

    # Create initial DataStruct with raw x and y
    d_raw = DataStruct(
        size(x, 2), size(x, 1), size(y, 2),
        zeros(1,1), zeros(1,1), zeros(1,1), zeros(1,1),
        0, 0, 0, 0,
        zeros(1,1), zeros(1,1), zeros(1,1), zeros(1,1),
        x, y, polynomial_order, data, file_name, test_ratio, x_features, y_features
    )

    # Return split and processed DataStruct
    return _split_data(d_raw)
end

function plot_data(d::DataStruct, i::Int, j::Int)
    x = d.Xtrain[:, i]
    y = d.Xtrain[:, j]
    labels = vec(d.Ytrain)

    scatter(x, y;
        group=labels,
        legend=false,
        xlabel="X$i",
        ylabel="X$j",
        markershape=:star,
        palette=:tab10,
        title="X$i vs X$j grouped by label"
    )
end

function plot_corr_row(d::DataStruct, k::Int)
    x = d.data[:, 1:end-1]
    y = vec(d.Y)
    nf = size(x, 2)

    plots = Vector{Any}(undef, nf)
    for i in 1:nf
        if i == k
            plots[i] = histogram(x[:, i]; title="X$i", xlabel="X$i")
        else
            plots[i] = scatter(x[:, k], x[:, i];
                group=y,
                markershape=:circle,
                legend=false,
                xlabel="X$k",
                ylabel="X$i",
                palette=:tab10
            )
        end
    end

    plot(plots..., layout=(1, nf), size=(250*nf, 250))
end

function plot_corr_matrix(d::DataStruct)
    x = d.data[:, 1:end-1]
    y = vec(d.Y)
    nf = size(x, 2)

    plots = Vector{Any}(undef, nf * nf)
    for i in 1:nf
        for j in 1:nf
            idx = (i - 1) * nf + j
            if i == j
                plots[idx] = histogram(x[:, i]; xlabel="X$i", title="X$i")
            else
                plots[idx] = scatter(x[:, j], x[:, i];
                    group=y,
                    xlabel="X$j",
                    ylabel="X$i",
                    legend=false,
                    markersize=3,
                    markershape=:circle,
                    palette=:tab10)
            end
        end
    end

    plot(plots..., layout=(nf, nf), size=(250 * nf, 250 * nf), title="Features Correlation Matrix")
end

function update_hyperparameter!(d::DataStruct, h::Dict{String,Any})
    if h["type"] == "testRatio"
        d.test_ratio = h["value"]
        _split_data(d)
    elseif h["type"] == "polyGrade"
        d.polynomial_order = h["value"]
    else
        @warn "Unknown hyperparameter type: $(h["type"])"
    end
end

function _generate_exponents(n_features::Int, target_deg::Int)
    current = zeros(Int, n_features)
    return _generate_exponents_recursive(n_features, target_deg, 1, current, Int[])
end

function _generate_exponents_recursive(
    n_features::Int, remaining_deg::Int, current_idx::Int,
    current::Vector{Int}, acc::Vector{Vector{Int}}
)
    if current_idx == n_features
        current[current_idx] = remaining_deg
        push!(acc, copy(current))
    else
        for i in 0:remaining_deg
            current[current_idx] = i
            _generate_exponents_recursive(n_features, remaining_deg - i, current_idx + 1, current, acc)
        end
    end
    return acc
end

function _split_data(d::DataStruct)
    nD = size(d.data, 1)
    TP = d.test_ratio
    r = randperm(nD)
    ntest = round(Int, TP / 100 * nD)
    ntrain = nD - ntest

    Xtrain = d.X[r[1:ntrain], :]
    Xtest  = d.X[r[ntrain+1:end], :]
    Ytrain = d.Y[r[1:ntrain], :]
    Ytest  = d.Y[r[ntrain+1:end], :]

    # Speed cubed
    Xtrain[:, 4] .= Xtrain[:, 4] .^ 3
    Xtest[:, 4]  .= Xtest[:, 4]  .^ 3

    # Wind direction as cosine
    Xtrain[:, 2] .= cosd.(Xtrain[:, 2])
    Xtest[:, 2]  .= cosd.(Xtest[:, 2])

    # Normalize X
    muX = mean(Xtrain, dims=1)
    sigmaX = std(Xtrain, dims=1)
    Xtrain_norm = (Xtrain .- muX) ./ sigmaX
    Xtest_norm  = (Xtest  .- muX) ./ sigmaX

    # Normalize Y
    muY = mean(Ytrain, dims=1)
    sigmaY = std(Ytrain, dims=1)
    Ytrain_norm = (Ytrain .- muY) ./ sigmaY
    Ytest_norm  = (Ytest  .- muY) ./ sigmaY

    return DataStruct(
        size(Xtrain_norm, 2), size(Xtrain_norm, 1), size(Ytrain_norm, 2),
        Xtrain_norm, Ytrain_norm, Xtest_norm, Ytest_norm,
        ntest, 0, 0, 0,
        muX, sigmaX, muY, sigmaY,
        d.X, d.Y, d.polynomial_order, d.data, d.file_name, d.test_ratio, d.x_features, d.y_features
    )
end

end