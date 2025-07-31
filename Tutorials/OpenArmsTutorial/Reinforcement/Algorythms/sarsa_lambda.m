function [w, rmse] = sarsa_lambda(env, policyFunction, getActiveTiles, params)
    % Unpack SARSA(λ) and tile coding parameters
    nFeatures = params.n_features;
    alpha = params.alpha;
    gamma = params.gamma;
    lambda = params.lambda;
    epsilon = params.epsilon;
    nEpisodes = params.n_episodes;

    gifFreq = params.gifFreq;

    % Initialize weights
    w = zeros(nFeatures, 1);

    % Weight history for gif
    weight_history = zeros(nFeatures, nEpisodes/gifFreq);
    
    % RMSE reference
    %w_ref_struct = load('w_reference.mat');
    %w_ref = w_ref_struct.weights;
    %rmse = zeros(nEpisodes, 1);

    for ep = 1:nEpisodes
        % Reset environment
        state = env.reset();
        [a, epsilon] = policyFunction(state, w, epsilon, params);

        % Initialize eligibility trace
        e = zeros(size(w));
        done = false;
        it = 0;

        % RMSE
        %tile_visit_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
        %tile_squared_errors = containers.Map('KeyType', 'char', 'ValueType', 'double');
        %step_count = 0;
        
        while ~done
            % Take step
            [next_state, reward, done] = env.step(state, a);

            % Choose next action
            [ap, epsilon] = policyFunction(next_state, w, epsilon, params);

            % Get active features
            idx = getActiveTiles(state, a, params);
            idx_p = getActiveTiles(next_state, ap, params);

            % Estimate Q-values
            Q = sum(w(idx));
            Qp = sum(w(idx_p));

            % RMSE:
            %{
                Q_ref = sum(w_ref(idx));
                step_count = step_count + 1;
                
                % Use tile index set as string key
                key = mat2str(idx);
                
                % Accumulate visitation count
                if isKey(tile_visit_counts, key)
                    tile_visit_counts(key) = tile_visit_counts(key) + 1;
                    tile_squared_errors(key) = tile_squared_errors(key) + (Q - Q_ref)^2;
                else
                    tile_visit_counts(key) = 1;
                    tile_squared_errors(key) = (Q - Q_ref)^2;
                end
            %}
            % TD error
            if done
                delta = reward - Q;
            else
                delta = reward + gamma * Qp - Q;
            end

            % Update eligibility trace (replacing)
            e(idx) = 1;

            % Update weights
            w = w + alpha * delta * e;

            % Decay traces
            e = gamma * lambda * e;

            % Transition
            state = next_state;
            a = ap;
            it = it + 1;

            if mod(ep, gifFreq) == 0
                weight_history(:, ep/gifFreq) = w;
            end
        end
        % Print progress
        disp(['Episode ', int2str(ep), ' completed in ', int2str(it), ' iterations'])

        
        % RMSE
        %{
        keys_ = keys(tile_visit_counts);
        weighted_error = 0;
        
        for i = 1:length(keys_)
            k = keys_{i};
            freq = tile_visit_counts(k) / step_count;  % d(s)
            err = tile_squared_errors(k);
            weighted_error = weighted_error + freq * err;
        end
        
        rmse(ep) = sqrt(weighted_error);
        %}
        
        %{
        % Store total return per episode
        initialState = [-0.5; 0]; maxSteps = 1000;
        greedyReturn = evaluate_greedy_episode(w, initialState, maxSteps, params, env, getActiveTiles);
        returns_per_episode(ep) = greedyReturn;
        %}
    end
    %{
    figure;
    plot(1:params.n_episodes, returns_per_episode, 'LineWidth', 1.5);
    xlabel('Episode');
    ylabel('Greedy Return from [-0.5, 0]');
    title('Performance of Learned Policy Over Time');
    grid on;
    drawnow;
    %}

    % Plot RMSE
    %plot_rmse_over_episodes(rmse)
    
    % Plot GIF of Qvalues
    %makeQgif(weight_history, params, 'GIF_Qvalues.gif')
    
    % Plot GIF of Policy
    %makePolicyGif(weight_history, params, @getActiveTiles, 'GIF_Policy.gif');
    
end
