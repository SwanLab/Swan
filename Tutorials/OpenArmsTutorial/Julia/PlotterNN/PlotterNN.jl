module PlotterNN

export PlotterNNStruct, plotBoundary, plotNetworkStatus, drawConfusionMat, drawSurfaceResults, image

using ..Network
using ..CostNN
using ..Data
using Plots
using LinearAlgebra
using StatsBase

mutable struct PlotterNNStruct
    data::Any
    network::Any
    costFunction::Any
    neuronsPerLayer::Vector{Int}
end

function PlotterNNStruct(s::Dict{String,Any})
    return PlotterNNStruct(
        s["data"],
        s["network"],
        s["costFunc"],
        s["network"].neuronsPerLayer
    )
end

function plotBoundary(obj::PlotterNNStruct, type::String) # !!! Cannot work without a method Data.buildModel()
    X = obj.data.Xtrain
    nF = size(X, 2)
    nPL = obj.neuronsPerLayer
    n_pts = 100
    graphzoom = 1

    x = createMesh(X, graphzoom, n_pts)
    h = computeHeights(obj, x[:,1], x[:,2], n_pts, nF)

    colorsc = ["r","g","b","c","m","y","k"]
    colorsc = reverse(colorsc[1:obj.data.nLabels])

    colorRGB = [ [1,0,0],
                 [0,1,0],
                 [0,0,1],
                 [0,1,1],
                 [1,0,1],
                 [1,1,0],
                 [0,0,0] ]
    colorRGB = reverse(colorRGB[1:obj.data.nLabels])

    if type == "contour"
        for i = 1:nPL[end]
            contour(x[:,1], x[:,2], h[:,:,i]', color=colorsc[i])
        end
    elseif type == "filled"
        for i = 1:nPL[end]
            heatmap(x[:,1], x[:,2], h[:,:,end+1-i]', alpha=0.3, color=colorsc[i])
        end
    elseif type == "filledS"
        I = mapslices(argmax, h; dims=3)
        rgb = zeros(n_pts, n_pts, 3)
        for i = 1:n_pts
            for j = 1:n_pts
                rgb[j,i,:] = colorRGB[I[i,j],:]
            end
        end
        heatmap(rgb[:,:,1], rgb[:,:,2], rgb[:,:,3], alpha=0.3)
    end

    title!("Contour 0")
    Data.plotdata(obj.data, 1, 2)
end

function createMesh(X, graphzoom, n_pts)
    extra_f1 = Base.mean(X[:,1]) * graphzoom
    extra_f2 = Base.mean(X[:,2]) * graphzoom
    x1 = LinRange(minimum(X[:,1])-extra_f1, maximum(X[:,1])+extra_f1, n_pts)
    x2 = LinRange(minimum(X[:,2])-extra_f2, maximum(X[:,2])+extra_f2, n_pts)
    x = hcat(x1, x2)
    return x
end

function plotNetworkStatus(obj::PlotterNNStruct)
    layer = obj.network.layer
    nPL = obj.neuronsPerLayer
    nLy = length(nPL)

    maxTH = maximum([maximum(abs.(l.W)) for l in layer])

    x_sep = 50
    y_sep = 30

    plt = plot(; xlim=(-20, nLy*x_sep-10), ylim=(-20, maximum(nPL)*y_sep+10),
                legend=false, xticks=[], yticks=[], aspect_ratio=1)

    neurons = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()

    for i in 1:nLy
        for j in 1:nPL[i]
            cx = (i-1) * x_sep
            cy = (j-1) * y_sep

            if i == 1
                color = :red
            elseif i == nLy
                color = RGB(0,0.45,0.74)
            else
                color = RGB(1,0.67,0.14)
            end

            draw_circle!(cx, cy, 10, color)
            neurons[(i,j)] = (cx, cy)
        end
    end

    # Connections
    for i in 1:nLy-1
        W = layer[i].W
        for j in 1:nPL[i]
            for k in 1:nPL[i+1]
                (bx, by) = neurons[(i,j)]
                (fx, fy) = neurons[(i+1,k)]

                wth = abs(W[j,k]) / maxTH
                lw = 3 * wth

                # Color mapping from red to blue based on weight
                linecolor = RGB(wth, 0, 1-wth)

                plot!([bx, fx], [by, fy], lw=lw, color=linecolor)
            end
        end
    end

    display(plt)
end

function drawConfusionMat(obj::PlotterNNStruct) # !!! Cannot work without getOutput()
    targets = obj.data.Ytest
    x = obj.data.Xtest
    outputs = obj.costFunction.getOutput(x)

    
    # Convert one-hot to label indices
    target_labels = findmax(targets, dims=2)[2][:]  # vector of true labels
    output_labels = findmax(outputs, dims=2)[2][:]  # vector of predicted labels

    # Build confusion matrix (Additional function needed because Julia does not have plotconfusion())
    cm = confusion_matrix(target_labels, output_labels)

    # Plot
    heatmap(cm, xlabel="Predicted", ylabel="True", title="Confusion Matrix")
end

function confusion_matrix(true_labels, pred_labels)
    classes = union(unique(true_labels), unique(pred_labels))
    n = length(classes)
    cm = zeros(Int, n, n)

    for (t,p) in zip(true_labels, pred_labels)
        cm[t,p] += 1
    end
    return cm
end

function drawSurfaceResults(obj::PlotterNNStruct) # Cannot work without getOutput()
    targets = obj.data.Ytest
    x = obj.data.Xtest
    outputs = obj.costFunction.getOutput(x)

    plotSurface(obj, targets', outputs')
end

function image(obj::PlotterNNStruct, row::Int) # Cannot work without getOutput()
    targets = obj.data.Ytest
    x = obj.data.Xtest
    outputs = obj.costFunction.getOutput(x)

    trg_vec = targets[row,:]
    out_vec = outputs[row,:]

    img_target = reshape(trg_vec, 28, 28)
    img_output = reshape(out_vec, 28, 28)

    # Show target image
    p1 = heatmap(img_target, color=:grays, title="Target Image", aspect_ratio=:equal)

    # Show output image
    p2 = heatmap(img_output, color=:grays, title="Output Image", aspect_ratio=:equal)

    display(p1)
    display(p2)
end

function computeHeights(obj::PlotterNNStruct, x1, x2, n_pts, nF) # !!! Cannot work without a method Data.buildModel()
    nPL = obj.neuronsPerLayer
    X_test = zeros(n_pts, nF, n_pts)
    h = zeros(n_pts * nPL[end], n_pts)
    h_3D = zeros(n_pts, n_pts, nPL[end])

    for i in 1:n_pts
        x2_aux = ones(n_pts) .* x2[i]
        xdata_test = hcat(x1, x2_aux)

        # Needs buildModel() in Data.jl (inexistant)
        #=
        tempData = DataStruct(X = xdata_test, polynomialOrder = obj.data.polyGrade)
        xful = Data.buildModel(tempData)

        X_test[:,:,i] = xful

        output = obj.costFunction.getOutput(X_test[:,:,i]) # Right now there's no getOutput method in CostNN. This might be outdated
        h[:,i] = reshape(output, n_pts * nPL[end])
        =#
    end
    #=
    for j in 1:nPL[end]
        h_3D[:,:,j] = h[(j-1)*n_pts+1 : j*n_pts, :]
    end

    return h_3D
    =#
end

function plotSurface(obj::PlotterNNStruct, target, output)
    nF = obj.data.nFeatures
    Ntest = size(target,2)

    result = zeros(nF, nF)

    for i = 1:Ntest
        j = findfirst(>=(0.5), output[:,i])
        k = findfirst(==(1), target[:,i])
        if !isnothing(j) && !isnothing(k)
            result[k,j] += 1
        end
    end

    quo = sum(diag(result))
    perc = (quo / Ntest) * 100

    surface(result, xlabel="Target", ylabel="Output", zlabel="Count", title="Output vs Expected Surface")
    println("Output value vs expected percentage is $perc %")
end

# Auxiliar function to draw circles for plotNetworkStatus()
function draw_circle!(x, y, r, color)
    θ = LinRange(0, 2π, 100)
    plot!(x .+ r*cos.(θ), y .+ r*sin.(θ), fill=(color,0.8), linecolor=color)
end

end