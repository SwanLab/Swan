function idx = getActiveTiles(state, a_idx, params)
    % Unpack
    stateRange = params.state_range;        % [n x 2]
    tilesPerDim = params.tiles_per_dim(:);  % [n x 1]
    nTilings = params.num_tilings;
    offsets = params.offsets;               % [n x nTilings]
    
    % Compute tile widths (n x 1)
    tileWidth = (stateRange(:,2) - stateRange(:,1)) ./ tilesPerDim;
    
    % Preallocate index output
    idx = zeros(nTilings, 1);
    %{
    % Precompute per-tiling base index offset
    if length(tilesPerDim) == 1
        dimProd = [1; tilesPerDim];
    else
        dimProd = [1; cumprod(tilesPerDim(1:end-1))];
    end
    %}
    dimProd = [1; cumprod(tilesPerDim(1:end-1))];
    tilePlaneSize = prod(tilesPerDim);
    baseOffset = (a_idx - 1) * nTilings * tilePlaneSize;
    
    for t = 1:nTilings
        % Shift state with offset
        shiftedState = state + offsets(:, t) .* tileWidth;
        
        % Tile coordinates in each dimension (clamped)
        coords = floor((shiftedState - stateRange(:,1)) ./ tileWidth) + 1;
        coords = max(1, min(coords, tilesPerDim));
        
        % Compute linear index in tile plane
        linearTile = 1 + (coords - 1)' * dimProd;
        
        % Compute full index into global weight vector
        idx(t) = baseOffset + (t - 1) * tilePlaneSize + linearTile;
    end
end
