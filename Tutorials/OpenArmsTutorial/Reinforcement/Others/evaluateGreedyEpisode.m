function episodeReturn = evaluateGreedyEpisode(weights, initialState, maxSteps, params, env, activeTiles)
    % Evaluates greedy policy from a fixed initial state
    % Returns the total reward (or number of steps if you prefer)

    nActions = params.n_actions;
    state = initialState;
    step = 1;
    totalReward = 0;

    while true
        % Compute Q-values
        qVals = zeros(nActions, 1);
        for a = 1:nActions
            idx = activeTiles.get(state, a);
            qVals(a) = sum(weights(idx));
        end
        [~, action] = max(qVals);

        % Take step
        [nextState, reward, done] = env.step(state, action);
        totalReward = totalReward + reward;

        state = nextState;
        step = step + 1;

        if done || step > maxSteps
            break;
        end
    end

    % Return either total reward or number of steps (negative reward = -1 per step)
    episodeReturn = totalReward;           % Total reward (typically negative)
    % episodeReturn = step - 1;            % Alternatively: number of steps (equivalent if reward = -1 per step)
end
