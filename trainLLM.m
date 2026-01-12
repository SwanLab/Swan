%% trainLLM.m
% Toy "LLM": character-level next-token (next-char) causal Transformer in MATLAB
% Robust on Linux/Windows/macOS (maps uint16 char codes -> ids)
% Avoids dlarray dimension labels to prevent pagemtimes label errors.
% Requires Deep Learning Toolbox: dlarray, dlgradient, dlfeval, adamupdate, pagemtimes.

clear; clc;

%% 1) Tiny dataset (edit freely)
lines = [
"i like pizza."
"i like sushi."
"you like pizza."
"we like food."
"transformers are cool."
"matlab can train models."
];
text = lower(join(lines, newline));   % one long string

%% 2) Build vocabulary + map chars -> ids (ROBUST)
txt = char(text);                     % char vector (includes newline char 10)

codes    = unique(uint16(txt));       % character codes present in corpus
V        = numel(codes);              % vocab size
codeToId = containers.Map(num2cell(codes), num2cell(1:V)); % code -> id
idToChar = char(codes);               % id -> char (index into this)

% Convert corpus to id sequence
ids = zeros(numel(txt),1,'int32');
for i = 1:numel(txt)
    ids(i) = int32(codeToId(uint16(txt(i))));
end

% Space id for padding
if isKey(codeToId, uint16(' '))
    spaceId = int32(codeToId(uint16(' ')));
else
    error("Your corpus has no space character; add one or change padding token.");
end

%% 3) Create training sequences (context length)
contextLen = 64;   % T
stride     = 1;

N = max(0, floor((numel(ids)-contextLen-1)/stride) + 1);
if N < 1
    error("Corpus too small for contextLen=%d. Reduce contextLen.", contextLen);
end

X = zeros(contextLen, N, 'int32'); % [T x N]
Y = zeros(contextLen, N, 'int32'); % [T x N]
for n = 1:N
    s = (n-1)*stride + 1;
    X(:,n) = ids(s : s+contextLen-1);
    Y(:,n) = ids(s+1 : s+contextLen);
end

%% 4) Model hyperparameters
dModel    = 64;
numHeads  = 4;
dFF       = 128;
numLayers = 2;
dropout   = 0.1;

if mod(dModel, numHeads) ~= 0
    error("dModel must be divisible by numHeads.");
end

%% 5) Initialize parameters (dlarray, unlabeled)
params = struct();

% token embedding and positional embedding
params.tokenEmbed = dlarray(0.02*randn(dModel, V, 'single'));              % [dModel x V]
params.posEmbed   = dlarray(0.02*randn(dModel, contextLen, 'single'));     % [dModel x T]

for L = 1:numLayers
    params.("Wq"+L) = dlarray(0.02*randn(dModel, dModel,'single'));
    params.("Wk"+L) = dlarray(0.02*randn(dModel, dModel,'single'));
    params.("Wv"+L) = dlarray(0.02*randn(dModel, dModel,'single'));
    params.("Wo"+L) = dlarray(0.02*randn(dModel, dModel,'single'));

    params.("W1"+L) = dlarray(0.02*randn(dFF, dModel,'single'));
    params.("b1"+L) = dlarray(zeros(dFF,1,'single'));
    params.("W2"+L) = dlarray(0.02*randn(dModel, dFF,'single'));
    params.("b2"+L) = dlarray(zeros(dModel,1,'single'));

    params.("ln1g"+L) = dlarray(ones(dModel,1,'single'));
    params.("ln1b"+L) = dlarray(zeros(dModel,1,'single'));
    params.("ln2g"+L) = dlarray(ones(dModel,1,'single'));
    params.("ln2b"+L) = dlarray(zeros(dModel,1,'single'));
end

params.Wout = dlarray(0.02*randn(V, dModel, 'single')); % [V x dModel]
params.bout = dlarray(zeros(V,1,'single'));

%% 6) Adam state init
trailingAvg = struct(); trailingAvgSq = struct();
names = fieldnames(params);
for i = 1:numel(names)
    trailingAvg.(names{i})   = zeros(size(params.(names{i})), 'like', params.(names{i}));
    trailingAvgSq.(names{i}) = zeros(size(params.(names{i})), 'like', params.(names{i}));
end

%% 7) Training loop
epochs    = 30;
batchSize = 32;
learnRate = 2e-3;
gradClip  = 1.0;

numBatches = max(1, ceil(N/batchSize));

for epoch = 1:epochs
    idx = randperm(N);
    epochLoss = 0;

    for b = 1:numBatches
        j = idx((b-1)*batchSize + 1 : min(b*batchSize, N));
        xBatch = X(:,j); % [T x B]
        yBatch = Y(:,j); % [T x B]

        % IMPORTANT: no dimension labels (avoids pagemtimes label error)
        xBatch = dlarray(single(xBatch));
        yBatch = dlarray(single(yBatch));

        [loss, grads] = dlfeval(@modelGradients, params, xBatch, yBatch, ...
            V, contextLen, dModel, numHeads, numLayers, dropout);

        % clip
        gnames = fieldnames(grads);
        for k = 1:numel(gnames)
            g = grads.(gnames{k});
            grads.(gnames{k}) = max(min(g, gradClip), -gradClip);
        end

        % Adam update
        beta1 = 0.9; beta2 = 0.999; eps = 1e-8;
        t = (epoch-1)*numBatches + b;

        for k = 1:numel(names)
            nme = names{k};
            [params.(nme), trailingAvg.(nme), trailingAvgSq.(nme)] = adamupdate( ...
                params.(nme), grads.(nme), trailingAvg.(nme), trailingAvgSq.(nme), ...
                t, learnRate, beta1, beta2, eps);
        end

        epochLoss = epochLoss + double(gather(extractdata(loss)));
    end

    fprintf("Epoch %d/%d - loss: %.4f\n", epoch, epochs, epochLoss/numBatches);
end

%% 8) Sampling / generation
prompt = "i like ";
genLen = 120;
temperature = 0.9;

out = generateText(params, prompt, genLen, temperature, ...
    codeToId, idToChar, spaceId, V, contextLen, dModel, numHeads, numLayers);

disp("=== Generated ===");
disp(out);

%% ========================= Helper functions =========================

function [loss, grads] = modelGradients(params, xTB, yTB, V, T, dModel, numHeads, numLayers, dropout)
    logits = forwardModel(params, xTB, V, T, dModel, numHeads, numLayers, dropout, true); % [V x T x B]

    % reshape for CE
    [~, Tdim, Bdim] = size(logits);
    logits2 = reshape(logits, V, Tdim*Bdim);
    labels  = reshape(yTB, 1, Tdim*Bdim);

    % Stable softmax-cross-entropy
    logits2 = logits2 - max(logits2,[],1);
    expz = exp(logits2);
    p = expz ./ sum(expz,1);

    lab = int32(gather(extractdata(labels)));
    cols = 1:numel(lab);
    idx = sub2ind([V, numel(lab)], double(lab), cols);
    ptrue = p(idx);

    loss = -mean(log(ptrue + 1e-12), 'all');
    grads = dlgradient(loss, params);
end

function logits = forwardModel(params, xTB, V, T, dModel, numHeads, numLayers, dropout, isTraining)
    % xTB: [T x B] token ids as dlarray(single), integer-valued
    xInt = int32(gather(extractdata(xTB))); % for embedding lookup on CPU
    [Tdim, Bdim] = size(xInt);

    tokE = gather(extractdata(params.tokenEmbed)); % [dModel x V]
    posE = gather(extractdata(params.posEmbed));   % [dModel x T]

    E = zeros(dModel, Tdim, Bdim, 'single');
    for b = 1:Bdim
        E(:,:,b) = tokE(:, xInt(:,b)) + posE(:, 1:Tdim);
    end

    h = dlarray(E); % [dModel x T x B], unlabeled

    % causal mask: -1e9 above diagonal
    m = triu(ones(Tdim,Tdim,'single'), 1);
    attnMask = -1e9 * m;

    for L = 1:numLayers
        h1 = layernormSimple(h, params.("ln1g"+L), params.("ln1b"+L));
        a  = mhaCausal(h1, params.("Wq"+L), params.("Wk"+L), params.("Wv"+L), params.("Wo"+L), ...
            numHeads, attnMask);

        if isTraining && dropout > 0
            a = dropoutSimple(a, dropout);
        end
        h = h + a;

        h2 = layernormSimple(h, params.("ln2g"+L), params.("ln2b"+L));
        f  = ffn(h2, params.("W1"+L), params.("b1"+L), params.("W2"+L), params.("b2"+L));

        if isTraining && dropout > 0
            f = dropoutSimple(f, dropout);
        end
        h = h + f;
    end

    logits = pagemtimes(params.Wout, h) + params.bout; % [V x T x B]
end

function y = mhaCausal(x, Wq, Wk, Wv, Wo, numHeads, attnMask)
    % x:  [dModel x T x B]
    % W*: [dModel x dModel]
    [dModel, T, B] = size(x);
    headDim = dModel / numHeads;
    scale = sqrt(headDim);

    Q = pagemtimes(Wq, x); % [dModel x T x B]
    K = pagemtimes(Wk, x);
    Vv= pagemtimes(Wv, x);

    % Split heads: [headDim x numHeads x T x B]
    Q = reshape(Q, headDim, numHeads, T, B);
    K = reshape(K, headDim, numHeads, T, B);
    Vv= reshape(Vv, headDim, numHeads, T, B);

    out = zeros(headDim, numHeads, T, B, 'like', x);

    % (Loop implementation for clarity/compatibility)
    for b = 1:B
        for h = 1:numHeads
            q = squeeze(Q(:,h,:,b)); % [headDim x T]
            k = squeeze(K(:,h,:,b)); % [headDim x T]
            v = squeeze(Vv(:,h,:,b));% [headDim x T]

            scores = (q' * k) / scale;  % [T x T]
            scores = scores + attnMask; % causal
            scores = scores - max(scores,[],2);
            A = exp(scores);
            A = A ./ sum(A,2);          % softmax rows

            o = v * A';                 % [headDim x T]
            out(:,h,:,b) = reshape(o, headDim, 1, T, 1);
        end
    end

    out = reshape(out, dModel, T, B);  % merge heads
    y = pagemtimes(Wo, out);
end

function y = ffn(x, W1, b1, W2, b2)
    y = pagemtimes(W1, x) + b1; % [dFF x T x B]
    y = relu(y);
    y = pagemtimes(W2, y) + b2; % [dModel x T x B]
end

function y = layernormSimple(x, g, b)
    % Normalize over channel dimension (dim 1)
    mu = mean(x,1);
    v  = mean((x - mu).^2,1);
    y  = (x - mu) ./ sqrt(v + 1e-5);
    y  = y .* g + b;
end

function y = dropoutSimple(x, p)
    mask = (rand(size(x),'like',x) > p);
    y = x .* mask / (1-p);
end

function out = generateText(params, prompt, genLen, temperature, codeToId, idToChar, spaceId, V, T, dModel, numHeads, numLayers)
    p = char(lower(prompt));
    ids = zeros(numel(p),1,'int32');
    for i = 1:numel(p)
        c = uint16(p(i));
        if isKey(codeToId, c)
            ids(i) = int32(codeToId(c));
        else
            ids(i) = spaceId;
        end
    end

    for t = 1:genLen
        x = int32(ones(T,1) * spaceId);
        take = min(T, numel(ids));
        x(end-take+1:end) = ids(end-take+1:end);

        xTB = dlarray(single(x)); % unlabeled
        logits = forwardModel(params, xTB, V, T, dModel, numHeads, numLayers, 0.0, false);
        logits = logits(:,end,1); % last position logits [V x 1]

        z = logits / temperature;
        z = z - max(z);
        pvec = exp(z);
        pvec = pvec / sum(pvec);

        pnp = gather(extractdata(pvec));
        nextId = sampleCategorical(pnp);
        ids(end+1,1) = int32(nextId);
    end

    charsOut = idToChar(double(ids));
    out = string(charsOut');
end

function k = sampleCategorical(p)
    c = cumsum(p(:));
    r = rand();
    k = find(c >= r, 1, 'first');
end
