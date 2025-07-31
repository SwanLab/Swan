function env = createMountainCarEnv(tile_params)
    % Fixed constants
    posRange = tile_params.pos_range;
    velRange = tile_params.vel_range;
    posMin = posRange(1);
    posMax = posRange(2);
    velMin = velRange(1);
    velMax = velRange(2);

    % Define reset function
    function state = reset()
        pos = posMin + (posMax - posMin) * rand();
        vel = velMin + (velMax - velMin) * rand();
        %pos = -0.5; vel = 0;
        state = [pos; vel];
    end

    % Define step function
    function [nextState, reward, done] = step(state, action)
        pos = state(1);
        vel = state(2);

        % Convert action: 1=left, 2=neutral, 3=right
        force = (action - 2) * 0.001;
        gravity = -0.0025 * cos(3 * pos);

        % Velocity update
        vel = vel + force + gravity;
        vel = min(max(vel, velMin), velMax);

        % Position update
        pos = pos + vel;
        pos = min(max(pos, posMin), posMax);

        if pos == posMin
            vel = 0;
        end

        nextState = [pos; vel];
        done = pos >= 0.5;
        reward = -1;
    end

    % Return environment struct
    env.reset = @reset;
    env.step = @step;
end
