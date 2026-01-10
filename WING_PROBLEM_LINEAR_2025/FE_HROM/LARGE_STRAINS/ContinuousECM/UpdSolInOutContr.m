function [xNEW,wNEW,ISNEGATIVE,ISOUT,VARCnew] =  UpdSolInOutContr(xNEW,wNEW,delta_q,DATA,VAR_SMOOTH_FE,POLYINFO,VARCnew)

if nargin == 0
    load('tmp1.mat')
end
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx

ISNEGATIVE = 0 ;
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
[ISousideifempty,~,POLYINFO_NEW]=     EvaluateBasisFunctionALL(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ;


% EXAMINE WHETHER THERE ARE TRANSITIONS OF POINTS BETWEEN ELEMENTS (AND IF THERE ARE "FORBIDDEN" TRANSITIONS)
 [VARCnew,ListForbiddenTransitions]  = ForbiddenTransitionsElements(POLYINFO,DATA,POLYINFO_NEW,VARCnew) ; 




NNEGATIVES  = length(find(wNEW<0)) ;


%Once we have the increment of the design variables, 
%we have to move the points within the domain, as well as introduce the pertinent variations of the weigths. 
%Then we have to check whether the constraints of the problem are observed. 
%Concerning the weights, numerical experience indicates that the existence of
%a few negative weights do not necessary lead to a solution with negative weights ---indeed, after some
% iterations, these weights may even dissapear, or in the next step, if they are small, they will be removed. 
% The user is to provide a threshold for then maximum number of negative weights that are allowed to exist.   


if all(wNEW>=0)  && ~isempty(ISousideifempty) && isempty(ListForbiddenTransitions)
    % disp('All points are admissible !!!!!!!')
    %  disp('')
    
else
    %dbstop('34')
    %  disp('Inadmissible  point ...')
    if any(wNEW<0)
        %   disp('Some of the weights are negative')
        %   INDLLL = find(wNEW<0);
        %   wNEG = wNEW(INDLLL) ;
        
        %   disp(['Negative weights = ',num2str(wNEG'/sum(wNEW)*100),' (%  volume)']) ;
        ISNEGATIVE = 1 ;
        disp(['Found ',num2str(length(find(wNEW<0))),' negative weights']);
        if NNEGATIVES > DATA.THRESHOLD_NUMBER_OF_NEGATIVE_WEIGHTS
            ISOUT = 1;
            disp(['Number of negative points above the imposed threshold'])
        end
        
    end
    
    if (isempty(ISousideifempty) &&  ISOUT == 0) || ~isempty(ListForbiddenTransitions)
        if ~isempty(POLYINFO_NEW.ListPointsOutside)
            disp(['Points: ', num2str(POLYINFO_NEW.ListPointsOutside'),' are outside   the domain ....'])
        end
        
        if DATA.FORCE_POINTS_TO_REMAIN_INSIDE ==1
            
            
       %     The other constraint is that the points must remain within the domain.
       %In this case, we cannot ignore this fact ---as done with the weights---
       %because the integrand basis functions are not defined outside the domain. 
       %In this scenario, we have found that the following improve overall convergence: 
       %it consists  returning  the points back to the domain, and ``freeze'' its position for the next iteration. 
            
            
            
            disp('Repeating the iteration by freezing the position of   points outside or in critical regions  ')
            ListPointsToFreeze = unique([ListForbiddenTransitions(:); POLYINFO_NEW.ListPointsOutside(:)]) ;
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
        else
            ISOUT = 1 ;
        end
    end
    
    
    
end




