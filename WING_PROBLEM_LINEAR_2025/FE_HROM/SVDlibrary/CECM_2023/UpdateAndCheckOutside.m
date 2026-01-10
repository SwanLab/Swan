function [xNEW,wNEW,ISOUT,VARCnew,POLYINFO_NEW] =  UpdateAndCheckOutside(xNEW,wNEW,delta_q,DATA,VAR_SMOOTH_FE,POLYINFO,VARCnew)

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




