function plotPolicyColorMap(weights, params, getActiveTiles)
    % Unpack parameters
    posRange = params.pos_range;
    velRange = params.vel_range;
    nActions = params.n_actions;
    Pixels1 = params.pixels1;
    Pixels2 = params.pixels2;

    % Grid points
    posVals = linspace(posRange(1), posRange(2), Pixels1);
    velVals = linspace(velRange(1), velRange(2), Pixels2);
    actionVals = zeros(Pixels2, Pixels1);

    for i = 1:Pixels1
        for j = 1:Pixels2
            pos = posVals(i);
            vel = velVals(j);
            state = [pos; vel];

            qVals = zeros(nActions, 1);
            for a = 1:nActions
                idx = getActiveTiles(state, a, params);
                qVals(a) = sum(weights(idx));
            end

            [~, bestA] = max(qVals);
            actionVals(j, i) = bestA - ceil(nActions / 2); % Maps e.g., 1→-1, 2→0, 3→1
        end
    end

    % Plot
    %fig = figure('visible','off');
    figure;
    imagesc(posVals, velVals, actionVals);
    set(gca, 'YDir', 'normal');
    xlabel('Position');
    ylabel('Velocity');
    title('Greedy Action per State');
    
    % Colormap setup (e.g., 3 discrete actions: left, neutral, right)
    if nActions == 3
        colormap([1 0 0; 1 1 0; 0 1 0]);  % Red, Yellow, Green
        clim([-1 1]);
        colorbar('Ticks', [-1 0 1], 'TickLabels', {'Left (-1)', 'Neutral (0)', 'Right (+1)'});
    else
        % For general action sizes, use colormap with N colors
        colormap(jet(nActions));
        clim([0, nActions - 1]);
        colorbar('Ticks', 0:(nActions - 1), 'TickLabels', string(0:(nActions - 1)));
    end
end
