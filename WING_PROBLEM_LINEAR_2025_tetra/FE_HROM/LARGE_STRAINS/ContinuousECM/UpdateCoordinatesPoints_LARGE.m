function [xNEW, wNEW, nF, ISNEGATIVE,POLYINFO,ISOUT  ] = UpdateCoordinatesPoints_LARGE(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO)

%dbstop('4')
if nargin == 0
    load('tmp2.mat')
end


%[PHIk_y,dPHIk_y,POLYINFO ]= EvaluateBasisFunctionAtX_approx(xNEW, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;
 [PHIk_y,dPHIk_y,POLYINFO]=     EvaluateBasisFunctionAtX_LARGE(xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO);


bNEW = PHIk_y'*wNEW ;
% Equation to be solved
% F =  (b-bNEW)
% Residual
Fk = b-bNEW ;
nF = norm(Fk)  ;
% Computation of the Jacobian matrix
D = [];
for idim = 1:length(dPHIk_y)
    DerPHIt_x = bsxfun(@times,dPHIk_y{idim},wNEW)' ;
    D = [D DerPHIt_x] ;    
end
% Therefore, the Jacobian matrix is finally
D = [D  PHIk_y'] ;

% Truncated SVD of Jacobian matrix
if DATALOC.USE_LEAST_NORM_SOLUTION == 1
    delta_q = lsqminnorm(D,Fk) ; 
else
    
    USE_SVDT_FILTER_BEFORE = 0 ;
    if USE_SVDT_FILTER_BEFORE == 1
        DATALOCSVD.RELATIVE_SVD = 1;
        TOL = 1e-10 ;
        [UU,SS,VV] = SVDT(D,TOL,DATALOCSVD) ;
        
        % Solution of the underdetermined system of equations (more unknowns than equation)
        if length(SS) == size(D,1)
            % No need to correct rank
            % -------------------------------------------------------
            delta_q = D\Fk ;    % Matlab produces a sparse solution (minimize  the l1 norm of delta_q)
        else
            disp(['Incomplete rank: number of equations = ',num2str(size(D,1)),'; rank = ',num2str(length(SS)),' (TOLrel =',num2str(TOL)  ,')']);
            UU = bsxfun(@times,UU',1./SS)'  ;
            delta_q = VV'\(UU'*Fk) ;
        end
    else
        delta_q = D\Fk ;
    end
    
end

%
% % Updated direction
% delta_q = D\Fk ;


%*************************************
% Leveraging the sparseness of delta_q
%*************************************
% ABANDONED !!! It offers only marginal benefits in terms of computational
% performance 
% See MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/06_HROM_homog2D/EnhancedCECM.mlx
% [npoints,ndim ]= size(xNEW);
% ndof = npoints*ndim ;  
% dx = delta_q(1:ndof) ; 
% dx = reshape(dx,ndim,[]) ; 
% dx_zeros = sum(dx,1) ;
% % IND_MOVINGe
% % These are the indexes of the 
% IND_MOVING = find(dx_zeros ~= 0 ) ; % Indexes moving points
% POLYINFO.PREVIOUS_STEP.IND_POINTS_CHANGE_POSITION = IND_MOVING ; % Indexes moving points
% POLYINFO.PREVIOUS_STEP.PHIk_y  = PHIk_y ; 
% POLYINFO.PREVIOUS_STEP.dPHIk_y  = dPHIk_y ; 
%  POLYINFO.PREVIOUS_STEP.ELEMENTS_CONTAINING_xNEW  = POLYINFO.ELEMENTS_CONTAINING_xNEW ; 

%nmoving = length(IND_MOVING)  
% dw = delta_q(ndof+1:end) ;
% IND_WEIGHTS_CHANGE = find(dw ~= 0 ) ; 


% See  MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/06_HROM_homog2D/EnhancedCECM.mlx

if DATALOC.FORCE_POINTS_TO_REMAIN_INSIDE == 0 
    [xNEW,wNEW,ISNEGATIVE,ISOUT] =  UpdateSolutionWP_NOTFORCED(xNEW,wNEW,delta_q,DATALOC,VAR_SMOOTH_FE,POLYINFO) ; 
else

[xNEW,wNEW,ISNEGATIVE,ISOUT] =  UpdateSolutionWP_FORCED(xNEW,wNEW,delta_q,DATALOC,VAR_SMOOTH_FE,POLYINFO) ; 
end