function DATA = ReshapeIterativeSnapshot(DATA,VAR,kiter)
%--------------------------------------------------------------------------
% DATA = ReshapeIterativeSnapshot(DATA, VAR, kiter)
%
% PURPOSE:
%   Trims the snapshot matrices stored in DATA.SNAP_ITER to retain only
%   the relevant iterations up to the final converged one (kiter).
%
% DESCRIPTION:
%   After the iterative solution process (e.g., Newton–Raphson) finishes
%   for a given time step, this function removes unused columns in the
%   DATA.SNAP_ITER structure. It keeps only the first (kiter+1) iterations,
%   corresponding to the actual iterations performed during convergence.
%
% INPUT:
%   DATA   : Structure containing the simulation state, including SNAP_ITER.
%   VAR    : Structure with current variables (used only for sizing).
%   kiter  : Index of the final iteration (0-based), such that we retain columns 1 to kiter+1.
%
% OUTPUT:
%   DATA   : Updated structure with truncated snapshot matrices.
%
% NOTE:
%   - This function assumes that SNAP_ITER was initialized and populated
%     during the iteration loop.
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp3.mat')
end

   if ~isempty(DATA.SNAP_ITER)
        % Initializing iterative snapshots (with the values of the previous time step )
        fff = fieldnames(DATA.SNAP_ITER) ;
        
        for iii = 1:length(fff)
            nrows = size(VAR.(fff{iii})) ;
            DATA.SNAP_ITER.(fff{iii}) = DATA.SNAP_ITER.(fff{iii})(:,1:kiter+1) ;
        end
    end