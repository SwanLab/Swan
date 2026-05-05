module Model

using Flux

export build_model

# =========================================================
# Deux architectures RÉELLEMENT différentes :
#
# :mlp — Flatten complet + couches denses
#
#   Feature maps (8, 8, 128, N)
#       → flatten → (8192, N)
#       → Dense(8192, 256, relu)
#       → Dropout
#       → Dense(256, 1)
#       → sigmoid
#
#   Principe : chaque pixel de chaque feature map devient
#   une entrée du Dense. Le réseau apprend "où" dans
#   l'image les features s'activent.
#   Paramètres : 8192*256 + 256 + 256*1 = ~2 millions
#   + Meilleur si la position spatiale de la fissure compte
#   - Sur-apprentissage facile sur petit dataset
#   - Lent (beaucoup de paramètres)
#
# :gap — Global Average Pooling + Dense direct
#
#   Feature maps (8, 8, 128, N)
#       → moyenne spatiale par feature map → (128, N)
#       → Dense(128, 1)
#       → sigmoid
#
#   Principe : chaque feature map est résumée par sa moyenne
#   sur toute l'image spatiale. Le réseau apprend "est-ce
#   que ce filtre s'active" sans mémoriser "où exactement".
#   Paramètres : 128*1 + 1 = 129
#   + Très peu de paramètres → bon sur petit dataset
#   + Invariant à la position de la fissure dans l'image
#   - Perd l'information spatiale
# =========================================================

function build_model(; classifier=:gap, dropout_rate=0.3)

    @assert classifier in (:mlp, :gap) "classifier doit être :mlp ou :gap"

    # Blocs convolutifs communs
    # Après 3x MaxPool(2,2) sur une image 64x64 :
    # feature maps de taille (8, 8, 128, N)
    backbone = [
        Conv((3,3), 1 => 32, relu; pad=1),
        GroupNorm(32, 8),
        MaxPool((2,2)),      # 64 → 32

        Conv((3,3), 32 => 64, relu; pad=1),
        GroupNorm(64, 8),
        MaxPool((2,2)),      # 32 → 16

        Conv((3,3), 64 => 128, relu; pad=1),
        GroupNorm(128, 8),
        MaxPool((2,2)),      # 16 → 8
                             # sortie : (8, 8, 128, N)
    ]

    if classifier == :mlp
        # Flatten complet : (8, 8, 128, N) → (8192, N)
        head = [
            Flux.flatten,                    # 8*8*128 = 8192
            Dense(8*8*128, 256, relu),
            Dropout(dropout_rate),
            Dense(256, 1),
            sigmoid,
        ]
        println("Classifieur : MLP  [Flatten(8192) → Dense(8192→256) → Dropout → Dense(256→1)]")
        println("             ~2M paramètres — recommandé avec > 2000 images")

    else  # :gap
        # Global Average Pooling : (8, 8, 128, N) → (128, N)
        # Moyenne spatiale sur chaque feature map indépendamment
        head = [
            GlobalMeanPool(),   # (8,8,128,N) → (1,1,128,N)
            Flux.flatten,       # (1,1,128,N) → (128,N)
            Dense(128, 1),
            sigmoid,
        ]
        println("Classifieur : GAP  [GlobalAvgPool → (128,N) → Dense(128→1)]")
        println("             ~129 paramètres — recommandé avec < 1000 images")
    end

    return Chain(backbone..., head...)
end

end