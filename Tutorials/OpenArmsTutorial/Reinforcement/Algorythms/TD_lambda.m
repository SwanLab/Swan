function w = TD_lambda(env, policyFunction, getActiveTiles, params, type)
    % type = 'sarsa' or 'qlearning'

    % Unpack parameters
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

        % Eligibility trace
        e = zeros(size(w));
        done = false;
        it = 0;

        while ~done
            % Take step
            [next_state, reward, done] = env.step(state, a);

            % Choose next action ε-greedily
            [ap, epsilon] = policyFunction(next_state, w, epsilon, params);

            % Feature indices
            idx = getActiveTiles(state, a, params);

            Q = sum(w(idx));

            switch type
                case 'sarsa'
                    idx_p = getActiveTiles(next_state, ap, params);
                    Qp = sum(w(idx_p));
                    if done
                        delta = reward - Q;
                    else
                        delta = reward + gamma * Qp - Q;
                    end

                    % Update eligibility trace (replacing)
                    e(idx) = 1;

                    % Weight update
                    w = w + alpha * delta * e;

                    % Decay trace
                    e = gamma * lambda * e;

                case 'qlearning'
                    Qp_all = zeros(1, nActions);
                    for a_i = 1:nActions
                        idx_p = getActiveTiles(next_state, a_i, params);
                        Qp_all(a_i) = sum(w(idx_p));
                    end
                    Qp_max = max(Qp_all);

                    if done
                        delta = reward - Q;
                    else
                        delta = reward + gamma * Qp_max - Q;
                    end

                    % Update eligibility trace (replacing)
                    e(idx) = 1;

                    % Weight update
                    w = w + alpha * delta * e;

                    % Check if chosen action was greedy
                    idx_ap = getActiveTiles(next_state, ap, params);
                    isGreedy = sum(w(idx_ap)) == Qp_max;

                    % Reset or decay trace
                    if done || ~isGreedy
                        e = zeros(size(e));
                    else
                        e = gamma * lambda * e;
                    end

                otherwise
                    error('Unknown type. Use ''sarsa'' or ''qlearning''.');
            end

            % Transition
            state = next_state;
            a = ap;
            it = it + 1;
        end

        % Optional printout
        disp(['Episode ', int2str(ep), ' completed in ', int2str(it), ' iterations']);
    end
end
