function [w, rmse] = sarsa_lambda(env, policyFunction, getActiveTiles, params)
    % Unpack SARSA(λ) and tile coding parameters
    nFeatures = params.n_features;
    alpha = params.alpha;
    gamma = params.gamma;
    lambda = params.lambda;
    epsilon = params.epsilon;
    nEpisodes = params.n_episodes;
    s0 = params.initial_state;

    %gifFreq = params.gifFreq;

    % Initialize weights
    w = zeros(nFeatures, 1);

    % Weight history for gif
    %weight_history = zeros(nFeatures, nEpisodes/gifFreq);
    
    returns_per_episode = zeros(nEpisodes, 1);

    sumQ = zeros(nEpisodes, 1);
    sumSteps = zeros(nEpisodes, 1);
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

            % Choose next action
            [ap, epsilon] = policyFunction(next_state, w, epsilon, params);

            % Get active features
            idx = getActiveTiles(state, a, params);
            idx_p = getActiveTiles(next_state, ap, params);

            % Estimate Q-values
            Q = sum(w(idx));
            Qp = sum(w(idx_p));
            
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
        %{
            if mod(ep, gifFreq) == 0
                weight_history(:, ep/gifFreq) = w;
            end
        %}
        end
        %{
        qMaxVals = evaluate_Q(w, params);
        sumQ(ep) = sum(sum(qMaxVals));
        stepsPerState = evaluate_steps(w, env, params);
        sumSteps(ep) = sum(sum(stepsPerState));
        %}
        % Print progress
        disp(['Episode ', int2str(ep), ' completed in ', int2str(it), ' iterations'])
    end
    %{
    fig1 = figure;
    plot(1:ep,-sumQ(1:ep), 'LineWidth', 1.5); hold on
    plot(1:ep, sumSteps(1:ep), 'LineWidth', 1.5)
    hold
    xlabel('Episode');
    title('Integral of Q over episodes');
    legend('Integral of Q values', 'Sum of R')
    grid on;
    drawnow;
    saveas(fig1, 'LinearQAndSteps.png');

    fig2 = figure;
    semilogy(1:ep,-sumQ(1:ep), 'LineWidth', 1.5); hold on
    semilogy(1:ep, sumSteps(1:ep), 'LineWidth', 1.5)
    hold
    xlabel('Episode');
    title('Integral of Q over episodes');
    legend('Integral of Q values', 'Sum of R')
    grid on;
    drawnow
    saveas(fig2, 'LogarythmicQAndSteps.png');
    %}
    % Plot GIF of Qvalues
    %makeQgif(weight_history, params, 'GIF_Qvalues.gif')
    
    % Plot GIF of Policy
    %makePolicyGif(weight_history, params, @getActiveTiles, 'GIF_Policy.gif');

    % Plot GIF of states
    %makeStepsGif(weight_history, env, params, 'GIF_Steps2.gif')
    
end
