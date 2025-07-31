function img = gif_plotQColor2D(weights, params, episode)
    posRange = params.pos_range;
    velRange = params.vel_range;
    nActions = params.n_actions;

    posVals = linspace(posRange(1), posRange(2), Pixels1);
    velVals = linspace(velRange(1), velRange(2), Pixels2);
    qMaxVals = zeros(Pixels2, Pixels1);  % rows = vel, cols = pos
    Pixels1 = params.pixels1;
    Pixels2 = params.pixels2;

    for i = 1:Pixels1
        for j = 1:Pixels2
            state = [posVals(i); velVals(j)];
            q_vals = zeros(nActions, 1);
            for a = 1:nActions
                idx = getActiveTiles(state, a, params);
                q_vals(a) = sum(weights(idx));
            end
            qMaxVals(j, i) = max(q_vals);  % (row j, col i)
        end
    end

    % Plot the heatmap and return the image
    fig = figure('visible','off');
    imagesc(posVals, velVals, -qMaxVals);
    axis xy;  % make y increase upward
    xlabel('Position'); ylabel('Velocity');
    title(sprintf('max Q(s,a) over state space, episode %d', episode));
    colorbar; colormap jet;
    set(gca, 'FontSize', 10);

    frame = getframe(fig);
    img = frame2im(frame);
    close(fig);
end
