"""
plot_sog.jl — Tracé de la vitesse (Speed Over Ground) en fonction du temps

Structure attendue :
    projet/
    ├── plot_sog.jl         ← ce script
    └── merge_csv/          ← dossier contenant les CSV
        ├── fichier1.csv
        └── ...

Lancement :
    julia plot_sog.jl

Le graphique sera sauvegardé à côté du script :
    projet/sog_plot.png
"""

using CSV
using DataFrames
using Dates
using Plots

# ─── Configuration ────────────────────────────────────────────────────────────

# ⚠️  Modifie ce nom pour changer le fichier à analyser
const FILENAME = "08_0130.csv"

# Dossier source (même emplacement que dans merge_csv.jl)
const SOURCE_FOLDER = joinpath(@__DIR__, "merge_csv/")

# Nom de la colonne vitesse (insensible à la casse)
const SOG_COLUMN = "SpeedOverGround (m/s)"

# Nom du fichier de sortie image
const OUTPUT_IMAGE = joinpath(@__DIR__, "sog_plot.png")

# ─── Fonctions ────────────────────────────────────────────────────────────────

function parse_date(s::AbstractString)
    s = strip(s)
    try
        return DateTime(s[1:19], dateformat"yyyy-mm-dd HH:MM:SS")
    catch
        try
            return DateTime(s, dateformat"yyyy-mm-dd H:MM:SS")
        catch
            return missing
        end
    end
end

function find_column(df::DataFrame, target::String)
    # Recherche insensible à la casse
    for col in names(df)
        if lowercase(strip(col)) == lowercase(strip(target))
            return col
        end
    end
    return nothing
end

# ─── Programme principal ──────────────────────────────────────────────────────

function main()
    filepath = joinpath(SOURCE_FOLDER, FILENAME)

    println("📂 Fichier source : $filepath")
    println()

    # ── Vérifications ─────────────────────────────────────────────────────────
    if !isfile(filepath)
        @error """
        Fichier introuvable : $filepath
        ➜ Vérifie la valeur de FILENAME en haut du script.
        """
        exit(1)
    end

    # ── Chargement ────────────────────────────────────────────────────────────
    println("📖 Chargement du CSV...")
    df = CSV.read(filepath, DataFrame;
        silencewarnings = true,
        missingstring   = ["", "NA", "N/A", "null", "NULL"],
    )
    println("   $(nrow(df)) lignes × $(ncol(df)) colonnes chargées")
    println("   Colonnes disponibles :")
    foreach(c -> println("     • $c"), names(df))
    println()

    # ── Colonne date (2ème colonne) ────────────────────────────────────────────
    date_col = names(df)[2]
    println("📅 Colonne date utilisée : \"$date_col\"")
    df[!, :_datetime_] = parse_date.(string.(df[!, date_col]))

    n_invalid_dates = sum(ismissing.(df[!, :_datetime_]))
    if n_invalid_dates > 0
        @warn "$n_invalid_dates ligne(s) avec date invalide — exclues du tracé"
    end

    # ── Colonne SOG ───────────────────────────────────────────────────────────
    sog_col = find_column(df, SOG_COLUMN)
    if sog_col === nothing
        @error """
        Colonne \"$SOG_COLUMN\" introuvable dans le fichier.
        Colonnes disponibles : $(join(names(df), ", "))
        ➜ Vérifie la valeur de SOG_COLUMN en haut du script.
        """
        exit(1)
    end
    println("🚤 Colonne vitesse utilisée : \"$sog_col\"")

    # ── Filtrage des lignes valides ────────────────────────────────────────────
    df_valid = filter(row -> !ismissing(row[:_datetime_]) && !ismissing(row[sog_col]), df)
    sort!(df_valid, :_datetime_)

    n_total  = nrow(df)
    n_valide = nrow(df_valid)
    println("   Lignes valides pour le tracé : $n_valide / $n_total")
    println()

    if n_valide == 0
        @error "Aucune donnée valide à tracer (dates et/ou SOG tous manquants)."
        exit(1)
    end

    # ── Statistiques rapides ──────────────────────────────────────────────────
    sog_vals = Float64.(df_valid[!, sog_col])
    println("📊 Statistiques SOG :")
    println("   Min    : $(round(minimum(sog_vals), digits=4))")
    println("   Max    : $(round(maximum(sog_vals), digits=4))")
    println("   Moyenne: $(round(mean(sog_vals),    digits=4))")
    println("   Écart-type : $(round(std(sog_vals), digits=6))")

    # Détection de valeurs suspectes (trop constantes)
    n_unique = length(unique(round.(sog_vals, digits=5)))
    println("   Valeurs uniques (arrondi 5 déc.) : $n_unique")
    if n_unique <= 3
        @warn "⚠️  Très peu de valeurs distinctes ($n_unique) — données potentiellement gelées !"
    end
    println()

    # ── Tracé ─────────────────────────────────────────────────────────────────
    println("🖼️  Tracé en cours...")

    times = df_valid[!, :_datetime_]
    gr()  # backend Plots

    p = scatter(
        times,
        sog_vals;
        title        = "Speed Over Ground — $(FILENAME)",
        xlabel       = "Temps",
        ylabel       = "SOG (nœuds)",
        legend       = false,
        markersize   = 2,
        markerstroke = stroke(0),
        color        = :steelblue,
        alpha        = 0.6,
        size         = (1200, 500),
        dpi          = 150,
        xrotation    = 30,
    )

    # Ligne de référence horizontale à la valeur médiane pour repérer les plateaux
    med_val = median(sog_vals)
    hline!(p, [med_val];
        color     = :red,
        linestyle = :dash,
        linewidth = 1,
        label     = "Médiane ($(round(med_val, digits=3)))",
        legend    = true,
    )

    savefig(p, OUTPUT_IMAGE)
    println("✅ Graphique sauvegardé : $OUTPUT_IMAGE")
    println()

    println("""
    ─────────────────────────────────────────
    ✅ Terminé !
       Fichier analysé : $FILENAME
       Image créée     : $OUTPUT_IMAGE
    ─────────────────────────────────────────
    Conseil : si tu vois une longue plage de points
    parfaitement alignés horizontalement, c'est
    probablement un problème d'acquisition (valeur gelée).
    """)
end

using Statistics
main()