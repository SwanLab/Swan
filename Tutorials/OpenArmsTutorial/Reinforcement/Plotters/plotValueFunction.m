function plotValueFunction(w, params)
    nActions = length(params.actions);
    nPoints = 200;  % points to plot along state_range
    X = linspace(params.state_range(1), params.state_range(2), nPoints);
    V = zeros(size(X));

    for i = 1:nPoints
        vals = zeros(nActions,1);
        for a_idx = 1:nActions
            % getActiveTiles should accept the whole params struct or relevant fields
            idx = getActiveTiles(X(i), a_idx, params);
            vals(a_idx) = sum(w(idx));
        end
        V(i) = max(vals);  % value of best action in state X(i)
    end

    figure;
    plot(X, V, 'LineWidth', 2);
    xlabel('Position x');
    ylabel('Estimated Value');
    title('Learned Value Function');
    grid on;
end
