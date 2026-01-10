function [G,b,sqW,DATA,y,yORIG,TOL,Gnorm,z] = SelectCand_and_expandG(G,W,DATA)
%--------------------------------------------------------------------------
% function [G, b, sqW, DATA, y, yORIG, TOL, Gnorm, z] = SelectCand_and_expandG(G, W, DATA)
%
% PURPOSE:
%   Expands the current matrix `G` (typically G = BasisFᵗ, already weighted 
%   by √W) by appending a normalized vector from the orthogonal complement 
%   of the current subspace. This procedure is relevant in greedy algorithms 
%   for empirical cubature or in adaptive reduced basis construction.
%
%   It also filters candidate points for integration according to user-defined 
%   tolerances and constraints.
%
% DESCRIPTION:
%   - Computes vector `a` orthogonal to the current subspace spanned by `G`.
%   - If the contribution of `a` (its "importance") is above a threshold,
%     it is appended as a new row in `G`.
%   - Computes the vector `b` as the integral of the basis functions.
%   - Filters out candidate integration points with low column norm in `G`.
%   - Optionally restricts to user-provided candidate points.
%
% INPUTS:
%   - G     : Matrix (nRows x M) of already √W-weighted basis functions (G = BasisFᵗ)
%   - W     : Vector of weights (length M)
%   - DATA  : Structure containing parameters:
%              - TOL                     : Tolerance for convergence
%              - TOLFilterCandidatePoints: Threshold for filtering small norm columns
%              - RemoveColumnsWithNegativeProjection: (unused here)
%              - IND_POINTS_CANDIDATES   : Subset of indices to restrict search
%
% OUTPUTS:
%   - G       : Possibly expanded matrix with new direction appended
%   - b       : Projection of √W onto the rows of G
%   - sqW     : Vector of √W values
%   - DATA    : Updated data structure
%   - y       : Filtered candidate indices
%   - yORIG   : Original list of y before potential filtering
%   - TOL     : Convergence tolerance
%   - Gnorm   : Column-wise 2-norms of G
%   - z       : Empty set for future selection of integration points
%
% REFERENCES:
%   - Continuous Empirical Cubature Method (CECM)
%   - Greedy sampling strategies in reduced-order modeling
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE
%--------------------------------------------------------------------------

% BasisF is already multiplied by sqrt(W)
% -----------------------------------------
% BasisF = G'
% Orthogonal complement of sqrt(W)
% See
% Expansion with the orthogonal complement
% ----------------------------------------
sqW = sqrt(W) ;
a = sqW - G'*(G*sqW ) ;  
V =sum(W) ;
aN = norm(a) ;
importanceA = aN^2/V ;
TOLLOC = 1e-10 ;
if importanceA > TOLLOC
    G = [ a'/aN;G] ;
end
% -------------------------------------------

b = G*sqW ;  % b Vector (exact integral)

Gnorm =sqrt(sum(G.*G,1)) ; % Norm of Columns of G
M = size(G,2) ;  % Number of FE points
DATA = DefaultField(DATA,'TOL',0) ; % Default tolerance for convergence
TOL = DATA.TOL ;
% INITIALIZATIONS
% ------------------------
z = [] ; % Set of integration points
y=1:M ;
DATA = DefaultField(DATA,'TOLFilterCandidatePoints',1e-6) ;
% GnormNOONE =sqrt(sum(G(1:end-1,:).*G(1:end-1,:),1)) ; % Norm of Columns of G(see old, wrong line above)
%GnormNOONE =sqrt(sum(G(2:end,:).*G(2:end,:),1)) ; % cHANGE, 23-mARCH-2022 
if DATA.TOLFilterCandidatePoints >0
    TOL_REMOVE = DATA.TOLFilterCandidatePoints*norm(b) ;
    %rmvpin = find(GnormNOONE(y)<TOL_REMOVE) ;
    rmvpin = find(Gnorm(y)<TOL_REMOVE) ;
    y(rmvpin) = [] ;
end
DATA = DefaultField(DATA,'RemoveColumnsWithNegativeProjection',0); % 4-Dec-2019

DATA = DefaultField(DATA,'IND_POINTS_CANDIDATES',[]) ;

if ~isempty(DATA.IND_POINTS_CANDIDATES)
    y = intersect(y,DATA.IND_POINTS_CANDIDATES) ;
end
yORIG = y ;

