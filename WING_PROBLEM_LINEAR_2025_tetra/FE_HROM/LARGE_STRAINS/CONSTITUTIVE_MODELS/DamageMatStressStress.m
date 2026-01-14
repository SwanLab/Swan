function C = DamageMatStressStress(S,nstrain)
% Compute matrix C = S*S^T, damage model, vectorized fashion
% JAHO, 22-Oct-2025, Balmes 185, Barcelona
if nargin == 0
    load('tmp2.mat')
    S = [S;S] ; 
end
% buildSeffDyads  Block-stacked outer products C_i = Seff_i * Seff_i^T
%
% INPUTS
%   S : (nstrain*ng) x 1 vector, stacked as [Seff_1; Seff_2; ...; Seff_ng]
%             where each Seff_i is nstrain x 1 in Voigt ordering.
%   nstrain   : number of Voigt stress components (3, 4, or 6).
%
% OUTPUT
%   C    : (nstrain*ng) x nstrain block-stacked matrix [C_1; C_2; ...; C_ng],
%             with each block C_i = Seff_i * Seff_i^T (size nstrain x nstrain).
%
% NOTES
% - Engineering-shear factors (the "2" in double contractions) do NOT appear
%   here: this is a pure dyadic product, not a contraction.
% - Works for 3D (nstrain=6), plane strain (nstrain=4), or plane stress (nstrain=3).
%
% Example:
%   % 3D case, two GPs:
%   % Seff_1 = [sxx1 syy1 szz1 sxy1 syz1 szx1].';
%   % Seff_2 = [sxx2 syy2 szz2 sxy2 syz2 szx2].';
%   S = [Seff_1; Seff_2];
%   C = buildSeffDyads(S, 6);
%   % C(1:6, :) == Seff_1 * Seff_1.'
%   % C(7:12, :) == Seff_2 * Seff_2.'

if ~isvector(S)
    error('SeffALL must be a column vector');
end
if ~ismember(nstrain,[3,4,6])
    error('ncomp must be 3, 4, or 6');
end
N = numel(S);
if mod(N, nstrain) ~= 0
    error('buildSeffDyads:SizeMismatch', ...
        'Length(S) must be a multiple of nstrain. Got %d and nstrain=%d.', N, nstrain);
end
ng = N / nstrain;

% Preallocate with same precision/class as input
C = zeros(N, nstrain, 'like', S);

% Build strided indices for each component across all Gauss points
base = 1:nstrain:N;  % starting indices of each GP block (rows)
% Column-by-column assembly (nstrain ≤ 6 so this is cheap and vectorized over ng)
for j = 1:nstrain
    % Take the j-th component from every GP: (ng x 1)
    Seff_j = S((j-1) + (1:nstrain:N));
    % Repeat each entry nstrain times to match rows in each block: (N x 1)
    Seff_j_rep = kron(Seff_j, ones(nstrain,1,'like',S));
    % Column j: elementwise product gives stacked blocks Seff_i * Seff_i(j)
    C(:, j) = S .* Seff_j_rep;
end
end
