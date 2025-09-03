function a = softmaxPolicy(state, w_a, params, getActiveTiles)
    nActions = params.n_actions;
    h = zeros(1, nActions);
    for a_i = 1:nActions
        idx = getActiveTiles(state, a_i, params);
        h(a_i) = sum(w_a(idx));
    end
    % Softmax probabilities
    exp_h = exp(h - max(h)); % subtract max for numerical stability
    probs = exp_h / sum(exp_h);
    a = find(rand <= cumsum(probs), 1, 'first');
end