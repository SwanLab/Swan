function w = qlearning_lambda(env, policyFunction, getActiveTiles, params)
    % Unpack Q(λ) and tile coding parameters
    nFeatures = params.n_features;
    alpha = params.alpha;
    gamma = params.gamma;
    lambda = params.lambda;
    epsilon = params.epsilon;
    nEpisodes = params.n_episodes;
    nActions = params.n_actions;

    % Initialize weights
    w = zeros(nFeatures, 1);

    for ep = 1:nEpisodes
        % Reset environment
        state = env.reset();
        [a, epsilon] = policyFunction(state, w, epsilon, params);

        % Initialize eligibility trace
        e = zeros(size(w));
        done = false;
        it = 0;
        while ~done
            % Take step
            [next_state, reward, done] = env.step(state, a);

            % Choose next action ε-greedily
            [ap, epsilon] = policyFunction(next_state, w, epsilon, params);

            % Get active features
            idx = getActiveTiles(state, a, params);
            Qp_all = zeros(1, nActions);
            
            for a_i = 1:nActions
                idx_p = getActiveTiles(next_state, a_i, params);
                Qp_all(a_i) = sum(w(idx_p));
            end

            Q = sum(w(idx));
            Qp_max = max(Qp_all);

            if done
                delta = reward - Q;
            else
                delta = reward + gamma * Qp_max - Q;
            end

            % Update eligibility trace (replacing)
            e(idx) = 1;

            % Update weights
            w = w + alpha * delta * e;

            % Determine if ap is greedy
            isGreedy = sum(w(getActiveTiles(next_state, ap, params))) == Qp_max;

            % Eligibility trace update
            if done || ~isGreedy
                e = zeros(size(e)); % reset traces if non-greedy action chosen
            else
                e = gamma * lambda * e; % decay traces
            end

            % Transition
            state = next_state;
            a = ap;
            it = it + 1;
        end

        % Print progress
        disp(['Episode ', int2str(ep), ' completed in ', int2str(it), ' iterations'])
    end
end
