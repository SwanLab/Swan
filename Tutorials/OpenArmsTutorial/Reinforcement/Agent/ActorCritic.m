classdef ActorCritic < handle

    properties (Access = private)
        e_c   % eligibility trace for critic
        e_a   % eligibility trace for actor
        w_c   % critic weights
        w_a   % actor weights
    end

    methods (Access = public)
        function obj = ActorCritic(params)
            obj.e_c = zeros(params.n_features/params.n_actions, 1); % n_features / n_actions
            obj.e_a = zeros(params.n_features, 1); % n_features
            obj.w_c = zeros(params.n_features/params.n_actions, 1);
            obj.w_a = zeros(params.n_features, 1);
        end

        function [w_c, w_a] = computeW(obj, idx_s, idx_s_next, idx_a, ...
                                       done, reward, gamma, ...
                                       alpha_c, alpha_a, lambda)
           
            delta = obj.computeDelta(done, reward, gamma, idx_s, idx_s_next);

            % --- Critic update ---
            obj.e_c(idx_s) = 1;
            obj.w_c = obj.w_c + alpha_c * delta * obj.e_c;
            obj.e_c = gamma * lambda * obj.e_c;

            % --- Actor update ---
            obj.e_a(idx_a) = 1;
            obj.w_a = obj.w_a + alpha_a * delta * obj.e_a;
            obj.e_a = gamma * lambda * obj.e_a;

            w_c = obj.w_c;
            w_a = obj.w_a;
        end
    end

    methods (Access = private)

        function delta = computeDelta(obj, done, reward, gamma, idx_s, idx_s_next)
            % Compute value estimates
            V     = sum(obj.w_c(idx_s));
            Vnext = sum(obj.w_c(idx_s_next));

            % TD error
            if done
                delta = reward - V;
            else
                delta = reward + gamma * Vnext - V;
            end
        end
    end
end
