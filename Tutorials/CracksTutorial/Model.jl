module Model

using Flux

export build_model

# =========================================================
# Deux architectures de classification disponibles :
#
# :mlp  — GlobalMeanPool → Dense(128,64) → Dropout → Dense(64,1)
#   + Plus de capacité d'expression (combine les features)
#   + Meilleur sur grand dataset (> 2000 images)
#   - Risque de sur-apprentissage sur petit dataset
#   Nb paramètres classifieur : 128*64 + 64*1 = 8 256
#
# :gap  — GlobalMeanPool → Dense(128,1)
#   + Minimal, chaque filtre vote directement
#   + Meilleur sur petit dataset (< 1000 images)
#   + Plus rapide à entraîner
#   - Moins de capacité à combiner les features entre elles
#   Nb paramètres classifieur : 128*1 = 128
#
# Les blocs convolutifs sont identiques dans les deux cas.
# =========================================================

function build_model(; classifier=:mlp, dropout_rate=0.3)

    @assert classifier in (:mlp, :gap) "classifier doit être :mlp ou :gap"

    # Blocs convolutifs communs aux deux architectures
    # GroupNorm au lieu de BatchNorm : stable quelle que soit
    # la taille du batch, indispensable sur petit dataset
    backbone = [
        # Bloc 1 — textures basiques (bords, gradients)
        Conv((3,3), 1 => 32, relu; pad=1),
        GroupNorm(32, 8),
        MaxPool((2,2)),            # 64 → 32

        # Bloc 2 — structures locales (lignes fines, fissures courtes)
        Conv((3,3), 32 => 64, relu; pad=1),
        GroupNorm(64, 8),
        MaxPool((2,2)),            # 32 → 16

        # Bloc 3 — patterns globaux (fissure traversante)
        Conv((3,3), 64 => 128, relu; pad=1),
        GroupNorm(128, 8),
        MaxPool((2,2)),            # 16 → 8

        # Réduit (8,8,128,N) → (128,N) par moyenne spatiale
        GlobalMeanPool(),
        Flux.flatten,
    ]

    # Classifieur selon le choix
    if classifier == :mlp
        head = [
            Dense(128, 64, relu),
            Dropout(dropout_rate),
            Dense(64, 1),
            sigmoid,
        ]
        println("Classifieur : MLP  [Dense(128→64) → Dropout($(dropout_rate)) → Dense(64→1)]")
    else  # :gap
        head = [
            Dense(128, 1),
            sigmoid,
        ]
        println("Classifieur : GAP  [Dense(128→1)] — minimal, adapté aux petits datasets")
    end

    return Chain(backbone..., head...)
end

end