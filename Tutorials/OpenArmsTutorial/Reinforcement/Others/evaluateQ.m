function Q = evaluateQ(weights, params, activeTiles)
    posRange = params.pos_range;
    velRange = params.vel_range;
    nActions = params.n_actions;
    Pixels1 = params.pixels1;
    Pixels2 = params.pixels2;

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
                idx = activeTiles.get(state, a, params);
                q_vals(a) = sum(weights(idx));
            end
            qMaxVals(j, i) = max(q_vals);  % (row j, col i)
        end
    end
    Q = qMaxVals;


