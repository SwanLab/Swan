function plotQSurface2D(weights, params, Pixels1, Pixels2, activeTiles)
    % Extract state ranges (assumes 2D state: [position; velocity])
    posRange = params.pos_range;
    velRange = params.vel_range;
    nActions = params.n_actions;

    % Discretize state space
    posVals = linspace(posRange(1), posRange(2), Pixels1);
    velVals = linspace(velRange(1), velRange(2), Pixels2);
    qMaxVals = zeros(Pixels2, Pixels1);

       % Loop over grid of states
    for i = 1:Pixels1
        for j = 1:Pixels2
            pos = posVals(i);
            vel = velVals(j);
            state = [pos; vel];
            q_vals = zeros(nActions, 1);
            for a = 1:nActions
                idx = activeTiles.get(state, a);
                q_vals(a) = sum(weights(idx));
            end
            qMaxVals(j, i) = max(q_vals);
        end
    end
    
    
    % Create meshgrid and surface plot
    figure;
    [P, V] = meshgrid(posVals, velVals);
    surf(P, V, -qMaxVals, 'EdgeColor', 'none');  % Use -Q for height
    xlabel('Position');
    ylabel('Velocity');
    zlabel('-max Q(s,a)');
    %title('Surface Plot of -max Q-value over State Space');
    title('Episode 4000')
    colorbar;
    colormap jet;
    view([-135 30]);
    hold off