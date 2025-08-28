function idx = getActiveTiles2D(state, a_idx, params)
    % state: [position; velocity]
    % a_idx: action index (1, 2, or 3)
    % stateRange: 2x2 matrix [pos_min pos_max; vel_min vel_max]
    pos_range = params.pos_range;
    vel_range = params.vel_range;
    stateRange = [pos_range; vel_range];
    tilesPerDim = params.tiles_per_dim;
    nTilings = params.num_tilings;
    offsets = params.offsets;
    
    tileWidth = (stateRange(:,2) - stateRange(:,1)) ./ tilesPerDim;
    idx = zeros(nTilings, 1);
    
    for t = 1:nTilings
        offset = offsets(:, t);
        
        shiftedState = state + offset .* tileWidth;
        tileCoords = floor((shiftedState - stateRange(:,1)) ./ tileWidth) + 1;
        tileCoords = max(1, min(tileCoords, tilesPerDim)); % clamp to tile limits
        
        %{
        % Normalization of state:
        state_norm = (state - stateRange(:,1)) ./ (stateRange(:,2) - stateRange(:,1));  % now in [0,1]
        shiftedState = state_norm + offset ./ tilesPerDim;  % normalized offsets
        tileCoords = floor(shiftedState .* tilesPerDim) + 1;
        tileCoords = max(1, min(tileCoords, tilesPerDim)); % clamp to tile limits
        %}
        % Linear index for tile in this tiling
        linearTile = sub2ind([tilesPerDim(2), tilesPerDim(1)], tileCoords(2), tileCoords(1));

        % Global feature index: store in a 1D vector
        idx(t) = (a_idx - 1) * nTilings * prod(tilesPerDim) + ...
                 (t - 1) * prod(tilesPerDim) + linearTile;
    end
end