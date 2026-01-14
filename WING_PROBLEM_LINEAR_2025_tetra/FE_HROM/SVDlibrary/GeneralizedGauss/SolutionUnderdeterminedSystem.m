function delta_q = SolutionUnderdeterminedSystem(DATALOC,D,Fk,w,m)
% Function created 19-Oct-2021. It returns the increment of the  
% variable containing both position and weights. q = [x;w]
% Before october 2021, it was simply delta_q = D\Fk ;
% Now it incorporates more sophisticated ways of solving the problem.
% METHOD = DATALOC.METHOD_SOLVE_UNDERDETERMINED_SYSTEM_EQUATIONS 
% If METHOD == 0. Default matlab option delta_q = D\Fk (which it turns out is the l1-lest norm solution)
%  If METHOD ==1. l1-least norm solution with positive constraints on the
%  weights, via linear programming 
%  METHOD ==2: Least-norm solution (L2 norm), via the pseudo-inverse 
if nargin == 0 
    load('tmp1.mat')
   % DATALOC.METHOD_SOLVE_UNDERDETERMINED_SYSTEM_EQUATIONS = 1 ; 
end
METHOD = DATALOC.METHOD_SOLVE_UNDERDETERMINED_SYSTEM_EQUATIONS ; 

if METHOD == 0 
    % It turns out that this produces a sparse solution 
    % -------------------------------------------------------
    DATALOCSVD.RELATIVE_SVD = 1;  
    TOL = 1e-10 ; 
    [UU,SS,VV] = SVDT(D,TOL,DATALOCSVD) ; 
 
    if length(SS) == size(D,1)
        % It turns out that this produces a sparse solution 
    % -------------------------------------------------------
        delta_q = D\Fk ;
        
    else
        disp(['Incomplete rank: number of equations = ',num2str(size(D,1)),'; rank = ',num2str(length(SS)),' (TOLrel =',num2str(TOL)  ,')']); 
        ELIMINATE_RANK_DEF =1;
        if ELIMINATE_RANK_DEF == 1
             
            UU = bsxfun(@times,UU',1./SS)'  ;
            delta_q = VV'\(UU'*Fk) ;
        else
            delta_q = D\Fk ;
        end
    end
    
elseif METHOD == 1
    
     % Least L1-norm solution with positive weights (via linear programming)
     ndim = size(D,2)/size(w,1) -1 ; 
     if ndim == 2
    delta_q =  LeastL1NormPosW(DATALOC,D,Fk,w,m)  ; 
     else
         error('Option not implemented')
     end
    
    
elseif  METHOD == 2 
    % -----------
    % 14-Oct-2021
    % -----------
   % error('Unreliable option. Better l1-norm')
    [UU,SS,VV] = svd(D,'econ') ; 
    SS = diag(SS) ; 
    disp('Using the SVD to update the solution')
    disp(['Total Number of equations = ',num2str(size(D,1))])
    disp(['Number of independent equations =  ',num2str(length(SS))]) ; 
    % D q = F  ; U S V^T q = F ;    q = V*S^{-1}*(U^T *F)
    VV = bsxfun(@times,VV',1./SS)' ; 
    delta_q = VV*(UU'*Fk) ; 
  %  delta_q_old_method = D\Fk ; 
    
   % diffMETHODS = norm(delta_q-delta_q_old_method)/norm(delta_q_old_method) ; 
else
    error('Option not implemented')
end
