function plotEpisodeTrajectory2D(weights, initialState, maxSteps, params, env, activeTiles)
    
    nActions = params.n_actions;
    trajectory = zeros(maxSteps, 2);
    state = initialState;
    step = 1;

    while true
        % Store current state
        trajectory(step, :) = state';

        % Choose greedy action
        qVals = zeros(nActions, 1);
        for a = 1:nActions
            idx = activeTiles.get(state, a);
            qVals(a) = sum(weights(idx));
        end
        [~, action] = max(qVals);

        % Take environment step

        [nextState, ~, done] = env.step(state, action);
        state = nextState;
        step = step + 1;

        if done || step > maxSteps
            break;
        end
    end

        % Trim excess
    trajectory = trajectory(1:step-1, :);
    disp(['Episode finished in ', num2str(step - 1), ' steps.'])

    % Plot trajectory
    figure;
    plot(trajectory(:,1), trajectory(:,2), 'b-', 'LineWidth', 2); hold on;
    plot(trajectory(1,1), trajectory(1,2), 'go', 'MarkerSize', 8, 'LineWidth', 2);  % Start
    plot(trajectory(end,1), trajectory(end,2), 'rx', 'MarkerSize', 8, 'LineWidth', 2);  % End
    xlabel('State Dimension 1');
    ylabel('State Dimension 2');
    title('Greedy Policy Episode Trajectory');
    axis([params.pos_range, params.vel_range]);
    grid on;
end






