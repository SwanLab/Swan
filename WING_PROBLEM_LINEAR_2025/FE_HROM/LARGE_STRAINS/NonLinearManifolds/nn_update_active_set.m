function [w_after,SOLUTION_FOUND] = nn_update_active_set(Uloc, b, w_before)
%==========================================================================
% nn_update_active_set
%
% PURPOSE
% -------
% Compute a **nonnegative, equality-feasible** update of weights w that stays
% as close as possible (in Euclidean norm) to an initial guess w_before:
%
%     minimize  ||w - w_before||_2
%     subject to Uloc' * w = b,   w >= 0 .
%
% This is used as a local redistribution step in MAW-ECM/CECM pipelines
% after tentatively removing one quadrature point.
%
% PROBLEM SIZE / SHAPES
% ---------------------
% Uloc : (n × mB) or (mB × n) ?  — In this implementation Uloc is used as
%        UF = Uloc(F,:) so Uloc must be **(n × mB)**, i.e. columns are
%        equality constraints and rows correspond to weights. Then UF' has
%        size (mB × nF).  (Equivalently: size(Uloc,2) = mB, length(w)=n.)
% b    : (mB × 1)
% w_*  : (n × 1)
%
% INPUTS
% ------
% Uloc      : local constraint matrix; equality is enforced as Uloc' * w = b.
% b         : right-hand side of equality constraints.
% w_before  : initial weights (serves both as warm start and proximity center).
%
% OUTPUTS
% -------
% w_after        : feasible nonnegative solution; empty if infeasible.
% SOLUTION_FOUND : logical flag indicating success.
%
% METHOD (ACTIVE-SET WITH LEAST-NORM CORRECTION)
% ----------------------------------------------
% 1) Initialize free set F = {i | w_before(i) > 0}. Tight set T = complement.
% 2) On F, solve the **least-change equality correction**:
%       Let w_F = w_before_F + x, with UF' * x = b - UF' * w_before_F .
%    We compute the minimum-norm x with:
%       x = lsqminnorm( UF', b - UF' * w_before_F ).
%    Set w(F) = w_before(F) + x and w(T) = 0.
% 3) If any components of w are negative, move those indices to T (i.e., set
%    F(neg)=false) and repeat Step 2.
% 4) Stop when either:
%      • w >= 0 and equality holds within tolerance (success), or
%      • the iteration counter exceeds n - mB, or UF' is too “wide”
%        (rank/shape check), or residual is too large (failure).
%
% NUMERICAL SAFEGUARDS
% --------------------
% • Residual check:  ||UF'*x - (b - UF'*w_before_F)||_2  < 1e-10 * ||b - UF'*w_before_F||_2 .
% • Rank/shape guard: if size(UF',1) > size(UF',2) (more equations than free
%   unknowns), or residual test fails, declare infeasible for this F.
% • On failure, SOLUTION_FOUND=false and w_after=[] (explicit signal to caller).
%
% INTERPRETATION
% --------------
% The update equals the **orthogonal projection of w_before onto the affine
% equality subspace**, followed by an **NNLS-style active-set pruning**:
% we clamp negatives, shrink F, and re-project until w>=0 or we give up.
%
% COMPLEXITY (PER ITERATION)
% --------------------------
% Solves one lsqminnorm on (mB × nF); practical cost O(min(mB,nF)^2 * max(mB,nF)).
% Worst-case up to (n - mB + 1) iterations as F shrinks.
%
% WHEN CAN IT FAIL?
% -----------------
% • Equality constraints incompatible with w>=0 near w_before.
% • Too few free variables (nF < mB) after clamping.
% • Ill-conditioning of UF' (try re-scaling columns/rows or enriching F).
%
% TUNING / TIPS
% -------------
% • A warmer start (w_before closer to feasible) reduces clamping iterations.
% • If failures are frequent, consider relaxing previous elimination choices,
%   re-scaling Uloc columns, or adding a tiny floor to w_before before the
%   first iteration to seed F.
%
% AUTHOR / DATE
% -------------
% JAHO — Joaquín Alberto Hernández Ortega
% 21-Sep-2025, Barcelona
% (Comments updated by ChatGPT-5 Thinking)
%==========================================================================

%==========================================================================
if nargin == 0
    load('tmpE2.mat')
end


n = length(w_before); % Number of entries of the new weight vector
F = (w_before > 0);      % Number of positive entries of w_before
w = w_before;  % We seek w such that UF'*w-b  = 0, with w as close to w_before as possible
mB = length(b) ;  %  number of contraints
it = 0;  % Number of iterations, before 1-Oct-2025 it was set to 1
SOLUTION_FOUND = true;
while it <= n - mB
    % Solve on free set: minimize ||w_F - w_before_F|| s.t. U_F' * w_F + U_T' * w_T = b
    % Here w_T (tight) are zero; enforce equality with least-norm on free vars
    UF = Uloc(F,:);  %UT = Uloc(~F,:);
    %    rhs = b - UT'*zeros(nnz(~F),1);
    % Least-norm solution on free vars that satisfies equality:
    % find w_F minimizing ||w_F - w_before_F|| with UF' * w_F = rhs
    % -> set x = w_F - w_before_F, Aeq = UF', beq = rhs - UF' * w_before_F
    beq = b - UF' * w_before(F);
    x = lsqminnorm(UF', beq);      % unconstrained step on free vars
    residual = norm( UF'*x - beq)  ;
    if residual >=1e-10*norm(beq) || size(UF',1) > size(UF',2)
        %   disp('No feasible solution found')
        SOLUTION_FOUND = false;
        break
    end
    
    
    
    w(F) = w_before(F) + x;        % update free weights
    w(~F) = 0;
    
    % If any negatives, clamp and shrink free set
    neg = (w < 0);
    if ~any(neg)
        break
    else
        it = it+1 ;
    end
    F(neg) = false;                % move negatives to tight set and repeat
end
 

%w_after = max(w,0);   
% See explanation in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/13_HOMOG.mlx
if any(w < 0)
    SOLUTION_FOUND = false;
    w_after = [] ; 
else
    w_after = w;  
end

end

