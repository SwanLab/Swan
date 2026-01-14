function [xNEW,wNEW,ISOUT,VARCnew,POLYINFO_NEW] =  UpdateAndCheckOutside2023(xNEW,wNEW,delta_q,DATA,VAR_SMOOTH_FE,POLYINFO,VARCnew)
%--------------------------------------------------------------------------
% function [xNEW, wNEW, ISOUT, VARCnew, POLYINFO_NEW] = UpdateAndCheckOutside2023( ...
%                         xNEW, wNEW, delta_q, DATA, VAR_SMOOTH_FE, POLYINFO, VARCnew)
%
% PURPOSE:
%   Performs an update of quadrature points and weights in the context of
%   empirical cubature (e.g. DECM/CECM) and checks:
%     · If updated points remain inside the domain.
%     · If weights remain positive or within allowed threshold of negativity.
%
%   If any point goes outside the mesh domain, the update is rejected and
%   the corresponding point is frozen in future iterations (only its weight
%   can be updated). The function also ensures that the number of design
%   variables remains sufficient to satisfy integrand constraints.
%
% INPUTS:
% -------
%   - xNEW        : [N × ndim] current spatial coordinates of quadrature points.
%   - wNEW        : [N × 1] current weights of quadrature points.
%   - delta_q     : [ndof × 1] update vector: [Δx; Δw].
%                   · First part updates positions (displacements).
%                   · Second part updates weights.
%   - DATA        : structure with control fields:
%                   · THRESHOLD_NUMBER_OF_NEGATIVE_WEIGHTS
%   - VAR_SMOOTH_FE : structure with FE mesh and interpolation info.
%   - POLYINFO    : structure with polynomial interpolation metadata.
%   - VARCnew     : structure with:
%                   · POINTSl: indices of points whose position can change.
%                   · POINTRp: indices of points with fixed position.
%                   · DOFl   : degrees of freedom associated to displacement.
%
% OUTPUTS:
% --------
%   - xNEW          : updated coordinates (if no constraint violated).
%   - wNEW          : updated weights (if no constraint violated).
%   - ISOUT         : flag = 1 if update is rejected.
%   - VARCnew       : updated VARCnew structure (points may be frozen).
%   - POLYINFO_NEW  : updated interpolation metadata.
%
% METHOD:
% -------
%   1. Apply the updates to positions (Δx) and weights (Δw).
%   2. Use EVALBASIS to check if all updated points remain inside mesh.
%   3. If any weight becomes negative:
%       · Allow small negatives up to threshold.
%       · Otherwise, set ISOUT=1.
%   4. If any point leaves the domain:
%       · Revert update.
%       · Freeze its position by removing it from POINTSl and adding to POINTRp.
%       · Re-check if design space has enough DoFs.
%   5. Return updated variables or revert to original state.
%
% REMARKS:
% --------
%   - This function is typically called within optimization routines (e.g.,
%     trust-region or Newton-like schemes) for DECM or CECM formulations.
%   - Negative weights are tolerated up to a user-defined threshold to
%     improve convergence robustness.
%   - Domain violations are strictly forbidden due to invalid shape function
%     evaluations outside FE mesh.
%
% DEPENDENCIES:
% -------------
%   · EVALBASIS (typically alias for EvaluateBasisFunctionDIRECTFIT2023)
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, April 2023
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx

%ISNEGATIVE = 0 ;
ISOUT = 0 ;
% Update DOFl
ndofLOC = length(VARCnew.DOFl) ;
ndim = size(xNEW,2) ;
dX_L = reshape(delta_q(1:ndofLOC),ndim,[]) ;
xOLD  = xNEW ;
wOLD = wNEW ;
POINTS_F = [VARCnew.POINTSl(:);VARCnew.POINTSRp(:)] ;
xNEW(VARCnew.POINTSl,:) = xNEW(VARCnew.POINTSl,:) + dX_L' ;
% Weights
wNEW(POINTS_F) = wNEW(POINTS_F)  + delta_q(ndofLOC+1:end) ;

DATA.OnlyCheckIfIsInside = 1;
[ISousideifempty,~,POLYINFO_NEW]=     EVALBASIS(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ;

NNEGATIVES  = length(find(wNEW<0)) ;

%Once we have the increment of the design variables,
%we have to move the points within the domain, as well as introduce the pertinent variations of the weigths.
%Then we have to check whether the constraints of the problem are observed.
%Concerning the weights, numerical experience indicates that the existence of
%a few negative weights do not necessary lead to a solution with negative weights ---indeed, after some
% iterations, these weights may even dissapear, or in the next step, if they are small, they will be removed.
% The user is to provide a threshold for then maximum number of negative weights that are allowed to exist.


if all(wNEW>=0)  && ~isempty(ISousideifempty) %&& isempty(ListForbiddenTransitions)
    
else
    if any(wNEW<0)
        disp(['Found ',num2str(length(find(wNEW<0))),' negative weights']);
        if NNEGATIVES > DATA.THRESHOLD_NUMBER_OF_NEGATIVE_WEIGHTS
            ISOUT = 1;
            disp(['Number of negative points above the imposed threshold']) ;
            POLYINFO_NEW =  POLYINFO ; 
        end
    end
    
    if (isempty(ISousideifempty) &&  ISOUT == 0)
        if ~isempty(POLYINFO_NEW.ListPointsOutside)
            disp(['Points: ', num2str(POLYINFO_NEW.ListPointsOutside'),' are outside   the domain ....'])
        end
        
        
        %     The other constraint is that the points must remain within the domain.
        %In this case, we cannot ignore this fact ---as done with the weights---
        %because the integrand basis functions are not defined outside the domain.
        %In this scenario, we have found that the following improve overall convergence:
        %it consists  returning  the points back to the domain, and ``freeze'' its position for the next iteration.
        
        disp('Repeating the iteration by freezing the position of   points outside or in critical regions  ')
        ListPointsToFreeze = unique( POLYINFO_NEW.ListPointsOutside(:)) ;
        POLYINFO_NEW.setElements(ListPointsToFreeze) = POLYINFO.setElements(ListPointsToFreeze) ; %  
        VARCnew.POINTSRp = [VARCnew.POINTSRp; ListPointsToFreeze] ; % We include such points in the list
        % of points for which only the weights can change (but not the position)
        % Next we repeat the iteration, which means that we ignore the
        % updates carried out before. Consequently, we have to remove
        % it also from the list of POINTSl
        for icandidate = 1:length(ListPointsToFreeze)
            III = find(VARCnew.POINTSl ==ListPointsToFreeze(icandidate)) ;
            VARCnew.POINTSl(III) = [] ;
        end
        % Condition that the number of columns of the Jacobian should
        % be equal or greater than the number of integrand functions
        ncols = (ndim+1)*length(VARCnew.POINTSl) + length(VARCnew.POINTSRp) ;
        if ncols < size(xNEW,1)
            ISOUT = 1;
            disp('Number of design variables below number of equations...Getting out')
        end
        xNEW = xOLD ;
        wNEW = wOLD ;
        
    end
    
    
    
end




