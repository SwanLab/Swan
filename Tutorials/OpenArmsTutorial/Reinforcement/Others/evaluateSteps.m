function stepsPerState = evaluateSteps(weights, env, params, activeTiles)
    % Grid resolution
    posRange = params.pos_range;
    velRange = params.vel_range;
    Pixels1 = params.pixels1;
    Pixels2 = params.pixels2;
    nActions = params.n_actions;
    maxSteps = params.max_steps;  % to prevent infinite loops

    % Prepare grid
    posVals = linspace(posRange(1), posRange(2), Pixels1);
    velVals = linspace(velRange(1), velRange(2), Pixels2);
    stepsPerState = zeros(Pixels2, Pixels1);  % (vel, pos)

    % For each (pos, vel) combination
    count = 0;
    for i = 1:Pixels1
        for j = 1:Pixels2
            % Set initial state
            s0 = [posVals(i); velVals(j)];
            s = env.reset(s0);  % resets environment to desired state
            done = false;
            stepCount = 0;

            while ~done && stepCount < maxSteps
                % Greedy policy
                q_vals = zeros(nActions, 1);
                for a = 1:nActions
                    idx = activeTiles.get(s, a, params);
                    q_vals(a) = sum(weights(idx));
                end
                [~, a] = max(q_vals);

                % Step
                [s, ~, done] = env.step(s, a);
                stepCount = stepCount + 1;
            end
            count = count + 1;
            stepsPerState(j, i) = stepCount;
            disp(['Run ', int2str(count), ' completed'])
        end
    end
end
