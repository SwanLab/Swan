classdef OneDEnv
    properties
        actions
        walk_step
        arbitrariness
    end

    methods
        function obj = OneDEnv(params)
            obj.actions = params.actions;
            obj.walk_step = params.walk_step;
            obj.arbitrariness = params.arbitrariness;
        end

        function state = reset(s0)
            if nargin == 2
                state = s0;
            else
                state = rand() * 2 - 1;  % Random in [-1, 1]
            end
        end

        function [next_state, reward, done] = step(obj, state, action_idx)
            action = obj.actions(action_idx);
            noise = action * normrnd(0, obj.arbitrariness);
            next_state = state + obj.walk_step * action + noise;
            next_state = max(-1, min(1, next_state));  % Clamp to [-1, 1]

            reward = -next_state^2;
            done = abs(next_state) < 0.02;
        end
    end
end
