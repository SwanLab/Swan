function [w_c, w_a] = ActorCriticFunc(env, policyFunction, getActiveTiles, params, agent)

    % Extract params
    alpha_c   = params.alpha_c;
    alpha_a   = params.alpha_a;
    gamma     = params.gamma;
    lambda    = params.lambda;
    nEpisodes = params.n_episodes;
    nFeatures = params.n_features;

    for ep = 1:nEpisodes
        state = env.reset();
        w_a_dummy = zeros(nFeatures, 1);
        
        try
            a = policyFunction(state, w_a_dummy, params, getActiveTiles); % find a better way to do this w.o. dummy
        catch
            error('Policy function incompatible with ActorCritic agent');
        end

        done = false;
        it = 0;

        while ~done
            % Step environment
            [next_state, reward, done] = env.step(state, a);

            % Indices for critic and actor
            idx_s = getActiveTiles(state, 1, params);   % critic features
            idx_s_next = getActiveTiles(next_state, 1, params);
            idx_a = getActiveTiles(state, a, params);   % actor features

            % Update weights using class method
            [w_c, w_a] = agent.computeW(idx_s, idx_s_next, idx_a, ...
                                        done, reward, gamma, ...
                                        alpha_c, alpha_a, lambda);

            % Next action
            a = policyFunction(next_state, w_a, params, getActiveTiles);

            % Transition
            state = next_state;
            it = it + 1;
        end

        disp(['Episode ', int2str(ep), ' completed in ', int2str(it), ' iterations']);
    end
end
