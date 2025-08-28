classdef Agent
    properties (Access = private)
        env
        policyFunction
        getActiveTiles
        params
        w   % weight vector
    end
    
    methods
        function obj = Agent(env, policyFunction, getActiveTiles, params)
            % Constructor
            obj.env = env;
            obj.policyFunction = policyFunction;
            obj.getActiveTiles = getActiveTiles;
            obj.params = params;
            obj.w = zeros(params.n_features, 1);
        end
        
        function w = SARSA(obj)
            % Extract params for readability
            alpha = obj.params.alpha;
            gamma = obj.params.gamma;
            lambda = obj.params.lambda;
            epsilon = obj.params.epsilon;
            nEpisodes = obj.params.n_episodes;

            for ep = 1:nEpisodes
                % Reset environment
                state = obj.env.reset();
                [a, epsilon] = obj.policyFunction(state, obj.w, epsilon, obj.params);

                % Eligibility trace
                e = zeros(size(obj.w));
                done = false;
                it = 0;

                while ~done
                    % Step
                    [next_state, reward, done] = obj.env.step(state, a);

                    % Next action
                    [ap, epsilon] = obj.policyFunction(next_state, obj.w, epsilon, obj.params);

                    % Current Q
                    idx = obj.getActiveTiles(state, a, obj.params);
                    Q = sum(obj.w(idx));

                    % Next Q
                    idx_p = obj.getActiveTiles(next_state, ap, obj.params);
                    Qp = sum(obj.w(idx_p));

                    if done
                        delta = reward - Q;
                    else
                        delta = reward + gamma * Qp - Q;
                    end

                    % Update eligibility trace
                    e(idx) = 1;

                    % Weight update
                    obj.w = obj.w + alpha * delta * e;

                    % Decay trace
                    e = gamma * lambda * e;

                    % Transition
                    state = next_state;
                    a = ap;
                    it = it + 1;
                end
                disp(['Episode ', int2str(ep), ' completed in ', int2str(it), ' iterations']);
            end
            w = obj.w;
        end
        
        function w = Q_learning(obj)
            % Extract params for readability
            alpha = obj.params.alpha;
            gamma = obj.params.gamma;
            lambda = obj.params.lambda;
            epsilon = obj.params.epsilon;
            nEpisodes = obj.params.n_episodes;
            nActions = obj.params.n_actions;

            for ep = 1:nEpisodes
                % Reset environment
                state = obj.env.reset();
                [a, epsilon] = obj.policyFunction(state, obj.w, epsilon, obj.params);

                % Eligibility trace
                e = zeros(size(obj.w));
                done = false;
                it = 0;

                while ~done
                    % Step
                    [next_state, reward, done] = obj.env.step(state, a);

                    % Next action (for exploration only)
                    [ap, epsilon] = obj.policyFunction(next_state, obj.w, epsilon, obj.params);

                    % Current Q
                    idx = obj.getActiveTiles(state, a, obj.params);
                    Q = sum(obj.w(idx));

                    % Max Q over all actions
                    Qp_all = zeros(1, nActions);
                    for a_i = 1:nActions
                        idx_p = obj.getActiveTiles(next_state, a_i, obj.params);
                        Qp_all(a_i) = sum(obj.w(idx_p));
                    end
                    Qp_max = max(Qp_all);

                    if done
                        delta = reward - Q;
                    else
                        delta = reward + gamma * Qp_max - Q;
                    end

                    % Update eligibility trace
                    e(idx) = 1;

                    % Weight update
                    obj.w = obj.w + alpha * delta * e;

                    % Greedy check for eligibility trace reset
                    idx_ap = obj.getActiveTiles(next_state, ap, obj.params);
                    isGreedy = sum(obj.w(idx_ap)) == Qp_max;
                    
                    if done || ~isGreedy
                        e = zeros(size(e));
                    else
                        e = gamma * lambda * e;
                    end

                    % Transition
                    state = next_state;
                    a = ap;
                    it = it + 1;
                end
                disp(['Episode ', int2str(ep), ' completed in ', int2str(it), ' iterations']);
            end
            w = obj.w;
        end
    end
end
