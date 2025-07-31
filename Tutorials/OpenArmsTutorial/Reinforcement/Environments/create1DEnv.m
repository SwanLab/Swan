function env = create1DEnv(params)
    % Creates a 1D environment where the agent moves in [-1, 1] and tries to reach x = 0
    env.state = [];  % initial state will be set on reset

    % Initialize function handles for the environment
    env.reset = @() reset();
    env.step = @(state, action) step(state, action, params);

    function state = reset()
        % Reset state to a random position in [-1, 1]
        env.state = rand() * 2 - 1;
        state = env.state;
    end

    function [next_state, reward, done] = step(state, action_idx, params)
        % Apply action and return next state, reward, and done flag
        actions = params.actions;
        walk_step = params.walk_step;
        arbitrariness = params.arbitrariness;
        action = actions(action_idx);

        % Transition
        next_state = state + walk_step * action + action * normrnd(0, arbitrariness);
        next_state = max(-1, min(1, next_state));  % Clamp to [-1, 1]
        env.state = next_state;

        % Reward: shaped to encourage reaching 0
        reward = -next_state^2;

        % Terminal condition (optional): within small threshold around 0
        done = abs(next_state) < 0.02;
    end
end
