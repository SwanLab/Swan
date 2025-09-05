function plotMultipleTrajectories2D(weights, initialStates, maxSteps, params, env, activeTiles)
    % initialStates: matrix of size (nTrajectories × 2)
    
    nActions = params.n_actions;
    nTrajectories = size(initialStates, 1);

    % Set up figure
    figure;
    hold on;
    colors = lines(nTrajectories); % Distinct colors for each trajectory

    for i = 1:nTrajectories
        state = initialStates(i, :)';
        trajectory = zeros(maxSteps, 2);
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

            % Step in the environment
            [nextState, ~, done] = env.step(state, action);
            state = nextState;
            step = step + 1;

            if done || step > maxSteps
                break;
            end
        end

        % Trim excess
        trajectory = trajectory(1:step-1, :);

        % Plot the trajectory
        plot(trajectory(:,1), trajectory(:,2), '-', ...
            'Color', colors(i,:), 'LineWidth', 2);
        plot(trajectory(1,1), trajectory(1,2), 'go', ...  % green circle
            'MarkerSize', 8, 'LineWidth', 2);
        plot(trajectory(end,1), trajectory(end,2), 'rx', ...  % red X
            'MarkerSize', 8, 'LineWidth', 2);


        disp(['Trajectory ', num2str(i), ' finished in ', ...
            num2str(step - 1), ' steps.']);
    end

    xlabel('Position');
    ylabel('Velocity');
    title('Greedy Policy Trajectories from Different Initial States');
    axis([params.pos_range, params.vel_range]);
    grid on;
end
