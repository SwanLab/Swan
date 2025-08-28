function img = gif_plotPolicyColor2D(weights, params, getActiveTiles, episode)
    % Unpack parameters
    posRange = params.pos_range;
    velRange = params.vel_range;
    nActions = params.n_actions;
    Pixels1 = params.pixels1;
    Pixels2 = params.pixels2;

    % Grid
    posVals = linspace(posRange(1), posRange(2), Pixels1);
    velVals = linspace(velRange(1), velRange(2), Pixels2);
    actionVals = zeros(Pixels2, Pixels1);  % For imagesc: [row, col] = [vel, pos]

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
            actionVals(j, i) = bestA;  % Index 1-based
        end
    end

    % Plot to invisible figure
    fig = figure('Visible', 'off');
    imagesc(posVals, velVals, actionVals);
    set(gca, 'YDir', 'normal');
    xlabel('Position');
    ylabel('Velocity');
    title(sprintf('Greedy Action per State, Episode %d', episode));

    % Set colormap
    if nActions == 3
        colormap([1 0 0; 1 1 0; 0 1 0]);  % Red, Yellow, Green
        clim([1 3]);  % since bestA is 1,2,3
        colorbar('Ticks', 1:3, 'TickLabels', {'Left', 'Neutral', 'Right'});
    else
        colormap(jet(nActions));
        clim([1 nActions]);
        colorbar('Ticks', 1:nActions);
    end

    % Get image and destroy figure
    frame = getframe(fig);
    img = frame2im(frame);
    close(fig);
end
