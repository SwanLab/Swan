module Model

using Flux

export build_model

function build_model()
    return Chain(
        # Bloc 1
        Conv((3,3), 1 => 32, relu; pad=1),
        BatchNorm(32),
        MaxPool((2,2)),             # 128 → 64

        # Bloc 2
        Conv((3,3), 32 => 64, relu; pad=1),
        BatchNorm(64),
        MaxPool((2,2)),             # 64 → 32

        # Bloc 3
        Conv((3,3), 64 => 128, relu; pad=1),
        BatchNorm(128),
        MaxPool((2,2)),             # 32 → 16

        Flux.flatten,               # 128*16*16 = 32768

        Dense(128*16*16, 256, relu),
        Dropout(0.5),               # régularisation importante
        Dense(256, 1),
        sigmoid
    )
end

end