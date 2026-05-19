module Model

using Flux
using Metalhead

export build_model

# =========================================================
# build_model — deux modes disponibles
#
# model_type = :custom
#   → CNN from scratch avec GroupNorm
#   → classifier = :gap ou :mlp (voir détails ci-dessous)
#   → img_size modifiable par l'utilisateur
#
# model_type = :resnet
#   → ResNet18 pré-entraîné sur ImageNet (Metalhead.jl)
#   → Phase 1 : feature extraction — backbone gelé
#   → Seule la tête Dense(512, 1) est entraînée
#   → img_size fixe 128×128 RGB
# =========================================================

function build_model(; model_type=:custom,
                       classifier=:gap,
                       dropout_rate=0.3,
                       img_size=(128,128))

    @assert model_type in (:custom, :resnet) "model_type doit être :custom ou :resnet"

    if model_type == :custom
        return build_custom(; classifier=classifier,
                              dropout_rate=dropout_rate,
                              img_size=img_size)
    else
        return build_resnet(; dropout_rate=dropout_rate)
    end
end

# =========================================================
# CNN CUSTOM — from scratch
# =========================================================

function build_custom(; classifier=:gap, dropout_rate=0.3, img_size=(128,128))

    @assert classifier in (:mlp, :gap) "classifier doit être :mlp ou :gap"

    # Calcul de la taille des feature maps après 3x MaxPool(2,2)
    # Utilisé uniquement pour :mlp
    h_out = img_size[1] ÷ 8
    w_out = img_size[2] ÷ 8
    flat_size = h_out * w_out * 128

    backbone = [
        Conv((3,3), 1 => 32, relu; pad=1),
        GroupNorm(32, 8),
        MaxPool((2,2)),

        Conv((3,3), 32 => 64, relu; pad=1),
        GroupNorm(64, 8),
        MaxPool((2,2)),

        Conv((3,3), 64 => 128, relu; pad=1),
        GroupNorm(128, 8),
        MaxPool((2,2)),
    ]

    if classifier == :gap
        head = [
            GlobalMeanPool(),
            Flux.flatten,
            Dense(128, 1),
            sigmoid,
        ]
        println("Modèle   : Custom CNN")
        println("Classif. : GAP [GlobalAvgPool → Dense(128→1)]")
        println("img_size : $(img_size[1])×$(img_size[2]) grayscale")

    else  # :mlp
        head = [
            Flux.flatten,
            Dense(flat_size, 256, relu),
            Dropout(dropout_rate),
            Dense(256, 1),
            sigmoid,
        ]
        println("Modèle   : Custom CNN")
        println("Classif. : MLP [Flatten($flat_size) → Dense(256→1)]")
        println("img_size : $(img_size[1])×$(img_size[2]) grayscale")
    end

    return Chain(backbone..., head...)
end

# =========================================================
# RESNET18 — transfer learning (phase 1 : feature extraction)
#
# Architecture :
#   ResNet18 backbone (poids ImageNet gelés)
#   → AdaptiveAvgPool → Flatten → (512,)
#   → Dense(512, 1) → sigmoid
#
# Pourquoi geler le backbone ?
#   Les poids pré-entraînés encodent des features génériques
#   (bords, textures, formes) utiles pour toute image.
#   On ne veut pas les détruire avec un LR trop élevé.
#   On entraîne seulement la tête pour adapter la sortie
#   aux fissures — beaucoup plus rapide et stable.
#
# Pourquoi Dense(512, 1) ?
#   ResNet18 produit 512 feature maps en sortie du backbone.
#   GlobalAvgPool réduit chaque map à 1 valeur → vecteur 512.
#   Dense(512, 1) + sigmoid = classificateur binaire.
# =========================================================

function build_resnet(; dropout_rate=0.3)

    println("Modèle   : ResNet18 pré-entraîné (Metalhead)")
    println("Mode     : Feature extraction — backbone gelé")
    println("img_size : 128×128 RGB (fixe)")

    # Charge ResNet18 avec poids ImageNet
    resnet = ResNet(18; pretrained=true)

    # Extrait le backbone sans la tête de classification ImageNet
    # ResNet18 de Metalhead : resnet.layers contient le backbone + classifier
    backbone = resnet.layers[1]   # feature extractor

    # Gèle tous les paramètres du backbone
    # freeze! empêche le gradient de se propager dans ces couches
    Flux.freeze!(backbone)

    # Nouvelle tête adaptée à la classification binaire
    head = Chain(
        GlobalMeanPool(),      # (4, 4, 512, N) → (1, 1, 512, N)
        Flux.flatten,          # → (512, N)
        Dense(512, 64, relu),
        Dropout(dropout_rate),
        Dense(64, 1),
        sigmoid
    )

    model = Chain(backbone, head)

    # Compte les paramètres entraînables (tête seulement)
    total     = sum(length(p) for p in Flux.params(model))
    trainable = sum(length(p) for p in Flux.params(head))
    println("Paramètres total     : $total")
    println("Paramètres entraînables (tête) : $trainable")

    return model
end

end