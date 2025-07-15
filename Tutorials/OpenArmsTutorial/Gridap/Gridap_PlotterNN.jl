module PlotterNN

export PlotterNNStruct, init_plotter_nn, plot_boundary, plot_network_status, draw_confusion_mat, draw_surface_results, plot_image

using ..Network
using ..CostNN
using ..Data
using Plots
using LinearAlgebra
using StatsBase

"""
    PlotterNNStruct

Mutable struct storing neural network, data, and plotting configuration.
"""
mutable struct PlotterNNStruct
    data::DataStruct              
    network::Net                     
    cost_function::CostNNStruct     
    neurons_per_layer::Vector{Int}
end

"""
    init_plotter_nn(params::Dict{String,Any})

Factory function to initialize a PlotterNNStruct from a parameter dictionary.
"""
function init_plotter_nn(params::Dict{String,Any})
    return PlotterNNStruct(
        params["data"],
        params["network"],
        params["costFunc"],
        params["network"].neuronsPerLayer
    )
end

"""
    plot_boundary(obj, type)

Plots the decision boundary of the network with specified type ("contour", "filled", "filledS").
"""
function plot_boundary(obj::PlotterNNStruct, type::String) # !!! Cannot work without a method Data.buildModel()
    X = obj.data.Xtrain
    nF = size(X, 2)
    nPL = obj.neurons_per_layer
    n_pts = 100
    graphzoom = 1

    x = create_mesh(X, graphzoom, n_pts)
    h = compute_heights(obj, x[:,1], x[:,2], n_pts, nF)

    colorsc = reverse(["r","g","b","c","m","y","k"][1:obj.data.nLabels])
    colorRGB = reverse([[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[0,0,0]][1:obj.data.nLabels])

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
        for i = 1:n_pts, j = 1:n_pts
            rgb[j,i,:] = colorRGB[I[i,j],:]
        end
        heatmap(rgb[:,:,1], rgb[:,:,2], rgb[:,:,3], alpha=0.3)
    end

    title!("Contour 0")
    Data.plotdata(obj.data, 1, 2)
end

"""
    create_mesh(X, graphzoom, n_pts)

Creates a mesh grid based on training data and zoom factor.
"""
function create_mesh(X, graphzoom, n_pts)
    extra_f1 = mean(X[:,1]) * graphzoom
    extra_f2 = mean(X[:,2]) * graphzoom
    x1 = LinRange(minimum(X[:,1])-extra_f1, maximum(X[:,1])+extra_f1, n_pts)
    x2 = LinRange(minimum(X[:,2])-extra_f2, maximum(X[:,2])+extra_f2, n_pts)
    return hcat(x1, x2)
end

"""
    plot_network_status(obj)

Plots a graphical representation of the neural network architecture and weights.
"""
function plot_network_status(obj::PlotterNNStruct)
    layer = obj.network.layer
    nPL = obj.neurons_per_layer
    nLy = length(nPL)

    maxTH = maximum([maximum(abs.(l.W)) for l in layer])

    x_sep, y_sep = 50, 30

    plt = plot(; xlim=(-20, nLy*x_sep-10), ylim=(-20, maximum(nPL)*y_sep+10),
                legend=false, xticks=[], yticks=[], aspect_ratio=1)

    neurons = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()

    for i in 1:nLy
        for j in 1:nPL[i]
            cx, cy = (i-1)*x_sep, (j-1)*y_sep
            color = i == 1 ? :red : i == nLy ? RGB(0,0.45,0.74) : RGB(1,0.67,0.14)
            draw_circle!(cx, cy, 10, color)
            neurons[(i,j)] = (cx, cy)
        end
    end

    for i in 1:nLy-1
        W = layer[i].W
        for j in 1:nPL[i], k in 1:nPL[i+1]
            (bx, by), (fx, fy) = neurons[(i,j)], neurons[(i+1,k)]
            wth = abs(W[j,k]) / maxTH
            lw = 3 * wth
            linecolor = RGB(wth, 0, 1-wth)
            plot!([bx, fx], [by, fy], lw=lw, color=linecolor)
        end
    end

    display(plt)
end

"""
    draw_confusion_mat(obj)

Draws the confusion matrix for test data predictions.
"""
function draw_confusion_mat(obj::PlotterNNStruct)
    targets = obj.data.Ytest
    x = obj.data.Xtest
    outputs = obj.cost_function.getOutput(x)

    target_labels = findmax(targets, dims=2)[2][:]
    output_labels = findmax(outputs, dims=2)[2][:]

    cm = confusion_matrix(target_labels, output_labels)
    heatmap(cm, xlabel="Predicted", ylabel="True", title="Confusion Matrix")
end

"""
    confusion_matrix(true_labels, pred_labels)

Generates a confusion matrix from true and predicted labels.
"""
function confusion_matrix(true_labels, pred_labels)
    classes = union(unique(true_labels), unique(pred_labels))
    n = length(classes)
    cm = zeros(Int, n, n)
    for (t,p) in zip(true_labels, pred_labels)
        cm[t,p] += 1
    end
    return cm
end




end