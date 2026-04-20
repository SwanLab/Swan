module PlotterNN

export PlotterNNStruct, init_plotter_nn, plot_cost, plot_predictions

using ..Network
using ..CostNN
using ..Data
using Plots

"""
    PlotterNNStruct

Struct immutable pour la visualisation des résultats d'entraînement.
"""
struct PlotterNNStruct
    data::DataStruct
    network::Net
    cost_function::CostNNStruct
end

"""
    init_plotter_nn(params)

Initialise un PlotterNNStruct depuis un dictionnaire de paramètres.
"""
function init_plotter_nn(params::Dict{String, Any})
    return PlotterNNStruct(
        params["data"],
        params["network"],
        params["costFunc"]
    )
end

"""
    plot_cost(cost_history, title)

Affiche la courbe de coût au fil des époques.
cost_history : vecteur de valeurs de coût par époque.
"""
function plot_cost(cost_history::Vector{Float64}; title::String="Fonction coût")
    plt = plot(
        1:length(cost_history),
        cost_history;
        xlabel    = "Époques",
        ylabel    = "Coût",
        title     = title,
        linewidth = 1.8,
        legend    = false,
        grid      = true,
        yscale    = :log10
    )
    display(plt)
    return plt
end

"""
    plot_predictions(obj, θ)

Affiche Ypred vs Ytest sur les données de test (dénormalisées).
Utile pour vérifier visuellement la qualité de la régression.
"""
function plot_predictions(obj::PlotterNNStruct, θ::Vector{Float64})
    Xtest = obj.data.Xtest
    Ytest = obj.data.Ytest

    Ypred = compute_output(obj.network, Xtest, θ)

    # Dénormalisation
    Ypred_denorm = vec(Ypred) .* vec(obj.data.sigmaY) .+ vec(obj.data.muY)
    Ytest_denorm = vec(Ytest) .* vec(obj.data.sigmaY) .+ vec(obj.data.muY)

    plt = scatter(
        Ytest_denorm, Ypred_denorm;
        xlabel    = "Valeurs réelles",
        ylabel    = "Valeurs prédites",
        title     = "Prédictions vs Réalité",
        legend    = false,
        markersize = 3,
        alpha     = 0.6,
        grid      = true
    )

    # Ligne diagonale idéale y = x
    lims = (min(minimum(Ytest_denorm), minimum(Ypred_denorm)),
            max(maximum(Ytest_denorm), maximum(Ypred_denorm)))
    plot!(plt, collect(lims), collect(lims);
          linecolor = :red, linestyle = :dash, linewidth = 1.5)

    display(plt)
    return plt
end

end