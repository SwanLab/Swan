module Data

export DataStruct, plotdata, plotCorrRow, plotCorrMatrix, updateHyperparameter!


using CSV, DelimitedFiles, Statistics, Random, LinearAlgebra #, Plots, StatsPlots

mutable struct DataStruct
    # Public properties
    nFeatures::Int
    nSamples::Int
    nLabels::Int
    Xtrain::Matrix{Float64}
    Ytrain::Matrix{Float64}
    Xtest::Matrix{Float64}
    Ytest::Matrix{Float64}
    Ntest::Int
    batchSize::Int
    Batch_nD::Int
    Batch_nB::Int
    muX::Vector{Float64}
    sigmaX::Vector{Float64}
    muY::Vector{Float64}
    sigmaY::Vector{Float64}
    # Private properties
    X::Matrix{Float64}
    Y::Matrix{Float64}
    polynomialOrder::Int
    data::Matrix{Float64}
    fileName::String
    testRatio::Float64
    xFeatures::Vector{Int}
    yFeatures::Vector{Int}
end

function DataStruct(cParams::Dict{String, Any})
    # Extract parameters
    fileName        = cParams["fileName"]
    testRatio       = cParams["testRatio"]
    polynomialOrder = cParams["polynomialOrder"]
    xFeatures       = cParams["xFeatures"]
    yFeatures       = cParams["yFeatures"]

    # Load data
    data = readdlm(fileName, ',', Float64)
    x = data[:, xFeatures]
    y = data[:, yFeatures]

    X = x
    Y = y
    
    # Create initial dummy values for everything else
    dummyX = zeros(1, 1)
    dummyY = zeros(1, 1)
    dummyV = zeros(1)
    d = DataStruct(
        0, 0, 0, dummyX, dummyY, dummyX, dummyY, 0, 0, 0, 0,
        dummyV, dummyV, dummyV, dummyV,
        X, Y, polynomialOrder, data, fileName, testRatio, xFeatures, yFeatures
    )
     # Finish initialization
    splitdata!(d)
    d.nFeatures = size(d.Xtrain, 2)
    d.nLabels   = size(d.Ytrain, 2)
    d.nSamples  = size(d.Xtrain, 1)

    return d
end

function plotdata(d::DataStruct, i::Int, j::Int)
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

function plotCorrRow(d::DataStruct, k::Int)
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

function plotCorrMatrix(d::DataStruct)
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

function updateHyperparameter!(d::DataStruct, h::Dict{String,Any})
    if h["type"] == "testRatio"
        d.testRatio = h["value"]
        splitdata!(d)
    elseif h["type"] == "polyGrade"
        d.polynomialOrder = h["value"]
        # buildModel has been removed from the matlab code (calls are commented and method is inexistant)
        #d.X = buildModel(d)    
    else
        @warn "Unknown hyperparameter type: $(h["type"])"
    end
end

#=
function buildModel(d::DataStruct)
    X = d.X
    k = d.polynomialOrder
    x1 = X[:, 1]
    x2 = X[:, 2]
    n = size(X, 1)
    nTerms = div((k + 1) * (k + 2), 2) - 1  # Number of terms in 2D poly (excluding degree 0)
    Xful = Matrix{Float64}(undef, n, nTerms)

    cont = 1
    for g in 1:k
        for a in 0:g
            Xful[:, cont] = (x2 .^ a) .* (x1 .^ (g - a))
            cont += 1
        end
    end
    d.X = Xful
    return Xful
end

=#
function generateExponents(nFeatures::Int, targetDeg::Int)
    currentExponents = zeros(Int, nFeatures)
    return generateExponentsRecursive(nFeatures, targetDeg, 1, currentExponents, Int[])
end

function generateExponentsRecursive(
    nFeatures::Int,
    remainingDeg::Int,
    currentIndex::Int,
    currentExponents::Vector{Int},
    acc::Vector{Vector{Int}}
)
    if currentIndex == nFeatures
        currentExponents[currentIndex] = remainingDeg
        push!(acc, copy(currentExponents))
    else
        for i in 0:remainingDeg
            currentExponents[currentIndex] = i
            generateExponentsRecursive(
                nFeatures,
                remainingDeg - i,
                currentIndex + 1,
                currentExponents,
                acc
            )
        end
    end
    return acc
end


function splitdata!(d::DataStruct)
    nD = size(d.data, 1)
    TP = d.testRatio
    r = randperm(nD)
    ntest = round(Int, TP / 100 * nD)
    ntrain = nD - ntest

    d.Xtrain = d.X[r[1:ntrain], :]
    d.Xtest  = d.X[r[ntrain+1:end], :]
    d.Ytrain = d.Y[r[1:ntrain], :]
    d.Ytest  = d.Y[r[ntrain+1:end], :]

    # Speed cubed
    d.Xtrain[:, 4] .= d.Xtrain[:, 4] .^ 3
    d.Xtest[:, 4]  .= d.Xtest[:, 4]  .^ 3

    # Wind direction as cosine
    d.Xtrain[:, 2] .= cosd.(d.Xtrain[:, 2]) # cosd already takes degrees
    d.Xtest[:, 2]  .= cosd.(d.Xtest[:, 2])

    # Normalize X
    d.muX = mean(d.Xtrain, dims=1)[:]
    d.sigmaX = std(d.Xtrain, dims=1)[:]
    d.Xtrain = (d.Xtrain .- d.muX') ./ d.sigmaX'
    d.Xtest  = (d.Xtest  .- d.muX') ./ d.sigmaX'

    # Normalize Y
    d.muY = mean(d.Ytrain, dims=1)[:]
    d.sigmaY = std(d.Ytrain, dims=1)[:]
    d.Ytrain = (d.Ytrain .- d.muY') ./ d.sigmaY'
    d.Ytest  = (d.Ytest  .- d.muY') ./ d.sigmaY'
    d.Ntest = ntest
end

end