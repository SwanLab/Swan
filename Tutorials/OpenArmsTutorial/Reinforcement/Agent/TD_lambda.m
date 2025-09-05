function w = TD_lambda(env, policyFunction, getActiveTiles, params, agent)
    % type = 'sarsa' or 'qlearning'

    % Unpack parameters
    nFeatures = params.n_features;
    alpha = params.alpha;
    gamma = params.gamma;
    lambda = params.lambda;
    epsilon = params.epsilon;
    nEpisodes = params.n_episodes;
    

    % Initialize weights
    w = zeros(nFeatures, 1);

    for ep = 1:nEpisodes
        % Reset environment
        state = env.reset();
        [a, epsilon] = policyFunction(state, w, epsilon, params, @getActiveTiles);

        % Eligibility trace
        done = false;
        it = 0;

        while ~done
            % Take step
            [next_state, reward, done] = env.step(state, a);

            % Choose next action ε-greedily
            [ap, epsilon] = policyFunction(next_state, w, epsilon, params, @getActiveTiles);

            % Feature indices
            idx = getActiveTiles(state, a, params);

            
            Q = sum(w(idx));
            idx_p = getActiveTiles(next_state, ap, params);           

            w = agent.computeW(done,reward,gamma,w,idx,idx_p,Q,alpha,lambda,next_state,params);


            % Transition
            state = next_state;
            a = ap;
            it = it + 1;
        end

        % Optional printout
        disp(['Episode ', int2str(ep), ' completed in ', int2str(it), ' iterations']);
    end
end
