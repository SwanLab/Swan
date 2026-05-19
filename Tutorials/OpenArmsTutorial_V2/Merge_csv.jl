"""
merge_csv.jl — Fusion et tri de fichiers CSV par date

Structure attendue :
    projet/
    ├── merge_csv.jl        ← ce script
    └── merge_csv/          ← dossier contenant tous les CSV
        ├── fichier1.csv
        ├── fichier2.csv
        └── ...

Lancement depuis VSCode ou terminal :
    julia merge_csv.jl

Le fichier de sortie sera créé à côté du script :
    projet/merged_output.csv
"""

using CSV
using DataFrames
using Dates

# ─── Configuration ────────────────────────────────────────────────────────────

# Dossier source (relatif au script)
const SOURCE_FOLDER = joinpath(@__DIR__, "merge_csv/")

# Nom du fichier de sortie (créé à côté du script)
const OUTPUT_FILE = joinpath(@__DIR__, "merged_output.csv")

# Format de date dans la 2ème colonne  ex: "2025-05-16 0:34:58"
const DATE_FORMAT = dateformat"yyyy-mm-dd H:MM:SS"

# ─── Fonctions ────────────────────────────────────────────────────────────────

function parse_date(s::AbstractString)
    s = strip(s)
    try
        # Avec microsecondes : 2025-05-16 01:34:11.348181
        return DateTime(s[1:19], dateformat"yyyy-mm-dd HH:MM:SS")
    catch
        try
            # Sans microsecondes : 2025-05-16 0:34:58
            return DateTime(s, dateformat"yyyy-mm-dd H:MM:SS")
        catch
            return missing
        end
    end
end

function load_csv(path::String)
    try
        df = CSV.read(path, DataFrame;
            silencewarnings = true,
            missingstring   = ["", "NA", "N/A", "null", "NULL"],
        )

        if nrow(df) == 0
            @warn "Fichier vide, ignoré : $(basename(path))"
            return nothing
        end

        if ncol(df) < 2
            @warn "Moins de 2 colonnes dans $(basename(path)), ignoré."
            return nothing
        end

        date_col = names(df)[2]
        df[!, :_datetime_] = parse_date.(string.(df[!, date_col]))

        n_invalid = sum(ismissing.(df[!, :_datetime_]))
        if n_invalid > 0
            @warn "$(basename(path)) : $n_invalid ligne(s) avec date invalide — conservées en fin de fichier"
        end

        return df

    catch e
        @error "Impossible de lire $(basename(path)) : $e"
        return nothing
    end
end

function glob_csv(dir::String)
    files = readdir(dir; join = true)
    return filter(f -> isfile(f) && lowercase(splitext(f)[2]) == ".csv", files)
end

# ─── Programme principal ──────────────────────────────────────────────────────

function main()
    # Dossier du script lui-même

    source_dir  = SOURCE_FOLDER
    output_file = OUTPUT_FILE

    println("📂 Dossier CSV source : $source_dir")
    println("💾 Fichier de sortie  : $output_file")
    println()

    # Vérifie que le dossier source existe
    if !isdir(source_dir)
        @error """
        Dossier '$SOURCE_FOLDER' introuvable à côté du script.
        Vérifie que la structure est bien :
            projet/
            ├── merge_csv.jl
            └── merge_csv/   ← doit exister
        """
        exit(1)
    end

    # Recherche des CSV (exclut le fichier de sortie si jamais il est là)
    csv_files = filter(
        f -> abspath(f) != abspath(output_file),
        glob_csv(source_dir)
    )

    if isempty(csv_files)
        @error "Aucun fichier CSV trouvé dans : $source_dir"
        exit(1)
    end

    println("📄 $(length(csv_files)) fichier(s) CSV trouvé(s) :")
    foreach(f -> println("   • ", basename(f)), sort(csv_files))
    println()

    # ── Chargement ────────────────────────────────────────────────────────────
    frames   = DataFrame[]
    ref_cols = nothing

    for path in sort(csv_files)
        df = load_csv(path)
        if df !== nothing
            push!(frames, df)
            println("✅ Chargé : $(basename(path)) — $(nrow(df)) lignes")

            # Mémorise les colonnes du premier fichier valide comme référence
            if ref_cols === nothing
                ref_cols = setdiff(names(df), ["_datetime_"])
            end
        end
    end

    if isempty(frames)
        @error "Aucun fichier CSV valide chargé."
        exit(1)
    end

    # Avertit si certains fichiers ont des colonnes différentes
    for (i, df) in enumerate(frames[2:end])
        cols = setdiff(names(df), ["_datetime_"])
        manquantes = setdiff(ref_cols, cols)
        suppl      = setdiff(cols, ref_cols)
        !isempty(manquantes) && @warn "Fichier $i : colonnes manquantes : $manquantes"
        !isempty(suppl)      && @warn "Fichier $i : colonnes supplémentaires : $suppl"
    end

    # ── Fusion ────────────────────────────────────────────────────────────────
    println("\n🔀 Fusion en cours...")
    merged = vcat(frames...; cols = :union)

    total_avant = nrow(merged)

    # ── Tri par date ──────────────────────────────────────────────────────────
    merged[!, :_sort_key_] = coalesce.(merged[!, :_datetime_], typemax(DateTime))
    sort!(merged, :_sort_key_)

    # ── Dédoublonnage (supprime les timestamps en double) ─────────────────────
    unique!(merged, :_datetime_)

    doublons = total_avant - nrow(merged)
    if doublons > 0
        println("🗑️  $doublons ligne(s) dupliquée(s) supprimée(s)")
    else
        println("✅ Aucun doublon de date détecté")
    end

    # ── Nettoyage des colonnes temporaires ────────────────────────────────────
    select!(merged, Not([:_datetime_, :_sort_key_]))

    # Remet les colonnes dans l'ordre du premier fichier
    ordered_cols = vcat(
        filter(c -> c in names(merged), ref_cols),
        filter(c -> c ∉ ref_cols, names(merged))
    )
    select!(merged, ordered_cols)

    # ── Écriture ──────────────────────────────────────────────────────────────
    println("\n💾 Écriture du fichier de sortie...")
    CSV.write(output_file, merged)

    println("""
    ─────────────────────────────────────────
    ✅ Terminé !
       Fichier créé  : $output_file
       Lignes totales: $(nrow(merged))
       Colonnes      : $(ncol(merged))
    ─────────────────────────────────────────
    """)
end

main()