function offsets = generateRandomOffsets(nDims, nTilings)
%GENERATERANDOMOFFSETS Generates normalized random offsets for tile coding.
%
% Inputs:
%   nDims     - number of dimensions
%   nTilings  - number of tilings
%
% Output:
%   offsets   - (nDims x nTilings) matrix with values in [0, 1)

    offsets = rand(nDims, nTilings); % Each in [0, 1)
end