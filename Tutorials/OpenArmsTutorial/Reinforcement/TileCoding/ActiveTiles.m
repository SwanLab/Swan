classdef ActiveTiles < handle

    properties (Access = private)
        stateRange    % [n x 2]
        tilesPerDim   % [n x 1]
        nTilings      % scalar
        offsets       % [n x nTilings]
        tileWidth     % [n x 1]
        dimProd       % cumulative products for indexing
        tilePlaneSize % number of tiles in one tiling
    end

    methods
        function obj = ActiveTiles(params)
            % Constructor: initialize from params
            obj.stateRange    = params.state_range;
            obj.tilesPerDim   = params.tiles_per_dim(:);
            obj.nTilings      = params.num_tilings;
            obj.offsets       = params.offsets;

            % Precompute useful quantities
            obj.tileWidth     = (obj.stateRange(:,2) - obj.stateRange(:,1)) ./ obj.tilesPerDim;
            obj.dimProd       = [1; cumprod(obj.tilesPerDim(1:end-1))];
            obj.tilePlaneSize = prod(obj.tilesPerDim);
        end

        function idx = get(obj, state, a_idx)
            % Return active tile indices for given state and action index

            idx = zeros(obj.nTilings, 1);
            baseOffset = (a_idx - 1) * obj.nTilings * obj.tilePlaneSize;

            for t = 1:obj.nTilings
                % Shift state with offset
                shiftedState = state + obj.offsets(:, t) .* obj.tileWidth;

                % Tile coordinates in each dimension (clamped)
                coords = floor((shiftedState - obj.stateRange(:,1)) ./ obj.tileWidth) + 1;
                coords = max(1, min(coords, obj.tilesPerDim));

                % Compute linear index in tile plane
                linearTile = 1 + (coords - 1)' * obj.dimProd;

                % Compute full index into global weight vector
                idx(t) = baseOffset + (t - 1) * obj.tilePlaneSize + linearTile;
            end
        end
    end
end
