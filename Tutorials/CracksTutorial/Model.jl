module Model

using Flux

export build_model

# =========================================================
# Pourquoi ces choix d'architecture ?
#
# PROBLÈME OBSERVÉ : Train Acc = 100%, Test scores = ~0.5
# C'est du sur-apprentissage : le modèle mémorise les 420
# images au lieu d'apprendre à détecter les fissures.
#
# CAUSE 1 — BatchNorm sur petit dataset
#   BatchNorm normalise par batch. Avec ~30 images/batch,
#   les stats (moyenne, variance) sont trop bruitées →
#   gradients instables → val_loss qui saute (0.44→6.44→1.44)
#   → le modèle finit par mémoriser plutôt qu'apprendre.
#   FIX : GroupNorm normalise DANS chaque image, indépendant
#   du batch. Stable dès 1 image.
#
# CAUSE 2 — Dense(8192, 256) avec 420 images
#   8192 * 256 = 2 097 152 paramètres dans UNE couche.
#   Avec 420 images, c'est 4994 paramètres par image.
#   Le réseau mémorise pixel par pixel au lieu de généraliser.
#   FIX : GlobalMeanPool réduit (8,8,128,N) → (128,N)
#   en faisant la moyenne spatiale. On garde "ce filtre
#   détecte-t-il une fissure" sans mémoriser "où exactement".
#   Paramètres : 128*64 = 8192 au lieu de 2 097 152.
# =========================================================

function build_model(; dropout_rate=0.3)

    return Chain(

        # Bloc 1 — textures basiques (bords, gradients de gris)
        Conv((3,3), 1 => 32, relu; pad=1),
        GroupNorm(32, 8),     # 8 groupes de 4 channels — stable sur petits batches
        MaxPool((2,2)),       # 64 → 32

        # Bloc 2 — structures locales (lignes fines, début de fissure)
        Conv((3,3), 32 => 64, relu; pad=1),
        GroupNorm(64, 8),     # 8 groupes de 8 channels
        MaxPool((2,2)),       # 32 → 16

        # Bloc 3 — patterns globaux (fissure traversante, réseau de fissures)
        Conv((3,3), 64 => 128, relu; pad=1),
        GroupNorm(128, 8),    # 8 groupes de 16 channels
        MaxPool((2,2)),       # 16 → 8

        # GlobalMeanPool : (8, 8, 128, N) → (1, 1, 128, N) → (128, N)
        # Réduit 8192 → 128 valeurs, évite la mémorisation spatiale
        GlobalMeanPool(),
        Flux.flatten,

        # Classificateur léger — peu de paramètres = meilleure généralisation
        Dense(128, 64, relu),
        Dropout(dropout_rate),   # 0.3 : régularise sans bloquer l'apprentissage
        Dense(64, 1),
        sigmoid
    )
end

end