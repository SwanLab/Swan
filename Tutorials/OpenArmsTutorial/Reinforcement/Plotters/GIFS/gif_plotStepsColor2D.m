function img = gif_plotStepsColor2D(weights, env, params, episode)
    posRange = params.pos_range;
    velRange = params.vel_range;
    Pixels1 = params.pixels1;
    Pixels2 = params.pixels2;
    posVals = linspace(posRange(1), posRange(2), Pixels1);
    velVals = linspace(velRange(1), velRange(2), Pixels2);
    
    numSteps = evaluate_steps(weights, env, params);
    % Plot the heatmap and return the image
    fig = figure('visible','off');
    imagesc(posVals, velVals, numSteps);
    axis xy;  % make y increase upward
    xlabel('Position'); ylabel('Velocity');
    title(sprintf('Number of steps over state space, episode %d', episode));
    colorbar; colormap jet;
    set(gca, 'FontSize', 10);

    frame = getframe(fig);
    img = frame2im(frame);
    close(fig);
end
