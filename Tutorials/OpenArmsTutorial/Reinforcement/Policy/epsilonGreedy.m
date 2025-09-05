function [a, epsilon] = epsilonGreedy(state, w, epsilon, params, activeTiles)
    % Extract from params
    nActions = params.n_actions;
    epsilonStrategy = params.epsilon_strategy;

    % Decide whether to explore or exploit
    if rand < epsilon
        % Exploration
        a = randi(nActions);
    else
        % Exploitation
        values = zeros(1, nActions);
        for a_i = 1:nActions
            idx = activeTiles.get(state, a_i);
            values(a_i) = sum(w(idx));
        end
        % Break ties randomly
        bestActions = find(values == max(values));
        a = bestActions(randi(length(bestActions)));
    end

    % Update epsilon
    switch lower(epsilonStrategy)
        case 'static'
            % No change

        case 'decay'
            decayRate = params.epsilon_decay; % e.g., 0.9999
            epsilon = epsilon * decayRate;

        case 'linear_decay'
            minEpsilon = params.min_epsilon;
            decayRate = params.epsilon_decay;
            epsilon = max(minEpsilon, epsilon - decayRate);

        otherwise
            error("Unknown epsilon_strategy: %s", epsilonStrategy);
    end
end