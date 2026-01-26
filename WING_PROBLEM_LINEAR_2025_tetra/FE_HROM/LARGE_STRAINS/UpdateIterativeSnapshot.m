function DATA = UpdateIterativeSnapshot(DATA,VAR,iterk)
%--------------------------------------------------------------------------
% DATA = UpdateIterativeSnapshot(DATA, VAR, iterk)
%
% PURPOSE:
%   Stores the current values of selected variables (e.g., displacements,
%   stresses) during each iteration of the nonlinear solver (e.g., Newton–Raphson),
%   enabling post-analysis of convergence behavior.
%
% DESCRIPTION:
%   This function is used to fill the DATA.SNAP_ITER structure, which holds
%   the snapshot matrices for each variable being tracked at every iteration
%   of the current load/time step.
%
% INPUT:
%   DATA   : Structure containing simulation parameters and the SNAP_ITER field
%            where iterative snapshots are stored.
%   VAR    : Structure with current state variables (e.g., DISP, PK2STRESS, etc.)
%   iterk  : Current nonlinear iteration index.
%
% OUTPUT:
%   DATA   : Updated structure with current values of VAR fields stored in
%            SNAP_ITER for iteration iterk.
%
% NOTE:
%   - The fieldnames of SNAP_ITER must match those in VAR.
%   - This function assumes that DATA.SNAP_ITER is initialized beforehand.
%--------------------------------------------------------------------------

   if ~isempty(DATA.SNAP_ITER)
        % Initializing iterative snapshots (with the values of the previous time step )
        fff = fieldnames(DATA.SNAP_ITER) ;
        
        for iii = 1:length(fff)
            nrows = size(VAR.(fff{iii})) ;
            DATA.SNAP_ITER.(fff{iii})(1:nrows,iterk) = VAR.(fff{iii}) ;
        end
    end