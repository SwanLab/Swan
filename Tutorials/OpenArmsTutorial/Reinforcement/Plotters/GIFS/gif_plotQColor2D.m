function img = gif_plotQColor2D(weights, params, episode)
    posRange = params.pos_range;
    velRange = params.vel_range;
    Pixels1 = params.pixels1;
    Pixels2 = params.pixels2;
    posVals = linspace(posRange(1), posRange(2), Pixels1);
    velVals = linspace(velRange(1), velRange(2), Pixels2);
    
    qMaxVals = evaluate_Q(weights, params);

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
