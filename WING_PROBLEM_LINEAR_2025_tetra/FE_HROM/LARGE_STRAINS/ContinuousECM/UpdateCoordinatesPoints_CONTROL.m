function [xNEW, wNEW, nF, ISNEGATIVE,POLYINFO,ISOUT,VARCnew  ] = UpdateCoordinatesPoints_CONTROL(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew)

% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx 
if nargin == 0
    load('tmp.mat')
end

% Evaluation of the function and the gradients at all points (THIS SHOULD BE REVISED)
% SIGNIFICANT COMPUTATIONAL SAVINGS MAY ACCRUE FROM OPTIMIZING THIS PART OF
% THE CODE (3-JAN-2022)
% --------------------------------------------------------
 [PHIk_y,dPHIk_y,POLYINFO]=     EvaluateBasisFunctionALL(xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO);
 
 
 
 % VAlue of the function at POINTSRp
%  if isempty(VARCnew.POINTSRp)
%      PHIk_y_Rp = 0 ; 
%  else
%      PHIk_y_Rp = VARCnew.PHIk_y(VARCnew.POINTSRp,:)  ; 
%  end
%  % Value of the function
%  % Value of the function at POINTSRpw
% PHIk_y_Rpw = VARCnew.PHIk_y(VARCnew.POINTSRpw,:)  ; 

% RESIDUAL 
% -------
% This remains exactly the same 
% ------------------------------
bNEW = PHIk_y'*wNEW ;
ndim = size(xNEW,2) ; 
% Equation to be solved
% F =  (b-bNEW)
% Residual
r = b-bNEW ;
nF = norm(r)  ;
% Computation of the Jacobian matrix (derivative with respect to X)
J = zeros(length(r),length(wNEW)*ndim);
for idim = 1:length(dPHIk_y)    
    J(:,idim:ndim:end) =  bsxfun(@times,dPHIk_y{idim},wNEW)' ;
end
% We are only interested in the components associated to POINTSl 
ndim = size(xNEW,2) ; 
DOFl = small2large(VARCnew.POINTSl,ndim) ; 
VARCnew.DOFl  = DOFl ; 
J = J(:,DOFl) ; 
% The other two blocks are formed by 
POINTS_F = [VARCnew.POINTSl(:);VARCnew.POINTSRp(:) ] ; 
J = [J,  PHIk_y(POINTS_F,:)'] ;
%  Block decomposition of the Jacobian 



% Truncated SVD of Jacobian matrix
if DATALOC.USE_LEAST_NORM_SOLUTION == 1
    delta_q = lsqminnorm(J,r) ; 
else
    
    USE_SVDT_FILTER_BEFORE = 1 ;
    if USE_SVDT_FILTER_BEFORE == 1
        DATALOCSVD.RELATIVE_SVD = 1;
        TOL = 1e-10 ;
        [UU,SS,VV] = SVDT(J,TOL,DATALOCSVD) ;
        
        % Solution of the underdetermined system of equations (more unknowns than equation)
        if length(SS) == size(J,1)
            % No need to correct rank
            % -------------------------------------------------------
            delta_q = J\r ;    % Matlab produces a sparse solution (minimize  the l1 norm of delta_q)
        else
            disp(['Incomplete rank: number of equations = ',num2str(size(J,1)),'; rank = ',num2str(length(SS)),' (TOLrel =',num2str(TOL)  ,')']);
            UU = bsxfun(@times,UU',1./SS)'  ;
            delta_q = VV'\(UU'*r) ;
        end
    else
        delta_q = J\r ;
    end
    
end

[xNEW,wNEW,ISNEGATIVE,ISOUT,VARCnew] =  UpdSolInOutContr(xNEW,wNEW,delta_q,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew) ;


%


% See  MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/06_HROM_homog2D/EnhancedCECM.mlx

% if DATALOC.FORCE_POINTS_TO_REMAIN_INSIDE == 0 
%     [xNEW,wNEW,ISNEGATIVE,ISOUT] =  UpdateSolutionWP_NOTFORCED(xNEW,wNEW,delta_q,DATALOC,VAR_SMOOTH_FE,POLYINFO) ; 
% else
% 
% [xNEW,wNEW,ISNEGATIVE,ISOUT] =  UpdateSolutionWP_FORCED(xNEW,wNEW,delta_q,DATALOC,VAR_SMOOTH_FE,POLYINFO) ; 
% end