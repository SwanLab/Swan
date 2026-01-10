function [z,w,ERROR_GLO,DATAOUT]= DiscreteEmpiricalCubatureMethod(G,W,DATA)
%--------------------------------------------------------------------------
% function [z, w, ERROR_GLO, DATAOUT] = DiscreteEmpiricalCubatureMethod(G, W, DATA)
%
% PURPOSE:
%   Implements a robust version of the Discrete Empirical Cubature Method (DECM)
%   for hyperreduction in reduced-order modeling. The function selects a sparse set
%   of integration points (rows of G) and computes associated positive weights w
%   such that the orthogonality constraints are satisfied up to a user-defined tolerance.
%
%   This version ensures convergence by allowing the candidate set to be extended
%   dynamically, including the complementary set if the initial one fails to yield
%   a valid cubature rule.
%
% INPUTS:
%   - G     : [m x n] Matrix of orthogonal basis functions evaluated at Gauss points.
%             Each row corresponds to a Gauss point; each column is a basis function.
%   - W     : [m x 1] Gauss quadrature weights for the full-order model.
%   - DATA  : Structure with configuration options:
%              * TOL (scalar): relative tolerance for residual ||r|| / ||b||.
%              * IND_POINTS_CANDIDATES: optional list of initial candidate indices.
%              * STORE_INFO_ITERATIONS: flag to store error and selection history.
%
% OUTPUTS:
%   - z          : Indices of selected Gauss points (subset of rows of G).
%   - w          : Associated positive cubature weights for the selected points.
%   - ERROR_GLO  : Vector of error values at each iteration.
%   - DATAOUT    : Struct with additional outputs:
%                   * HistoryPoints
%                   * HistoryWeights
%                   * Error
%                   * kiteraciones (number of iterations)
%
% METHOD OVERVIEW:
%   1. Define the right-hand side vector b = Gᵗ * W.
%   2. Iteratively select rows (points) that reduce the residual r = b - Gᵗw.
%   3. Update the weights via rank-one inverse updates of the Gram matrix.
%   4. Ensure weights remain positive — if a weight becomes negative, the point is removed.
%   5. If convergence stalls (low progress after N iterations), enlarge the candidate set
%      by adding the complement (i.e., previously unused points).
%   6. Stop when the residual is below the given tolerance or max iterations are reached.
%
% FEATURES:
%   - Guaranteed convergence by fallback to full candidate space.
%   - Efficient weight update using inverse Gram matrix techniques.
%   - Robust to poor initial point selection via adaptive fallback.
%   - Tracks and stores history of selected points, weights, and residual errors (optional).
%
% REMARKS:
%   - This implementation corresponds to the enhanced DECM strategy described in
%     Hernández et al. (2024), improving upon the original NNLS-based cubature rules.
%   - The weights `w` are scaled by sqrt(W) to be consistent with L²(W)-orthogonality.
%   - The internal residual and objective are normalized to improve numerical stability.
%
% REFERENCES:
%   - Hernández et al. (2024), "CECM: A continuous empirical cubature method with application
%     to the dimensional hyperreduction of parameterized finite element models",
%     *Comput. Methods Appl. Mech. Engrg.*, Vol. 418, Article 116552.
%
% SEE ALSO:
%   - GetBasisMatrixECM2
%   - ContinuousECMgen2023
%   - UpdateWeightsInverse
%   - MultiUpdateInverseHermitian
%   - SelectCand_and_expandG
%
% AUTHOR:
%   J.A. Hernández (UPC/CIMNE), robust version of DECM as of June 2022.
% Comments by ChatGPT4
%--------------------------------------------------------------------------



% Version created 23-FEB-2021 to cope with clustering
% The only difference with respect to
% ../SVDlibrary/ECM/EmpiricalCubatureMethod_orig.m
% is that this version never fails to converge (in theory). If the
% set of candidates does not contain any solution, and the algorithm fails
% to converge in a given number of iterations, then the set of candidates is enlarged
% with the complementary, initial set, and then the search goes on...
% Modified 6-nov-2021
% See EmpiricalCubatureMethodCOMPL_aux.mlx
% Last modification: 21-June-2022. Removing elements whose entrines in G
% are small should be done before the expansion of the matrix
% ----------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end

%---------------------------------------------------------
% Expand G matrix and select candidate points
% ---------------------------------------------
[G,b,sqW,DATA,y,yORIG,TOL,Gnorm,z] = SelectCand_and_expandG(G,W,DATA) ;

DATA = DefaultField(DATA,'STORE_INFO_ITERATIONS',0) ;
if DATA.STORE_INFO_ITERATIONS == 1
    mpointsMAX  = size(G,1) ;
    DATAOUT.HistoryPoints = cell(1,mpointsMAX) ;
    DATAOUT.HistoryWeights = cell(1,mpointsMAX) ;
    DATAOUT.Error = zeros(1,mpointsMAX) ;
end
% -----------------------
% INITIALIZATIONS
% -----------------------

alpha = [] ; % Vector of weights
mPOS = 0 ; % Number of nonzero weights
r = b ; % Residual vector
k = 1;  % Number of iterations
errorGLO = [] ; %(for storing error)
% Default number of points
m =  length(b)  ;

normB = norm(b) ;
nerror = norm(r)/normB  ;
H = [] ; % Inverse of (Gz'*Gz)
ERROR_GLO = [] ;
NPOINTS =[] ;
nerrorACTUAL = nerror;
y = y(:) ;
NITERACIONES = 50*m ;
iOLD = [] ;
yCOMPL = (1:length(W))';
yCOMPL(y) = [] ; % Complementary set of the candidate set (in case this set is not sufficiently rich)


NITERATIONS_NO_MORE_POINTS = 10 ; %    NUMBER OF ITERATIONS IN WHICH THE CODE IS ALLOWED TO "ITERATE" IN ORDER TO FIND
%  THE MINIMUM. In Nov-4th-2021, a modification was introduced in this
%  regard, see EmpiricalCubatureMethod_CANDcompl_aux.mlx
nerrorMIN = 1e10 ;
COMPLEMENTARY_ALREADY_ADDED = 0 ;

ITER_MINIMUM = 0;
MAX_NUMBER_POINTS = 0 ;
% while  nerrorACTUAL >TOL && mPOS <m   && ~isempty(y) && k<NITERACIONES
% ...Before Nov. 4th -2021, see EmpiricalCubatureMethod_CANDcompl_aux.mlx
while  nerror >TOL && mPOS <m     && k<NITERACIONES
    
    % TIMELOC_k = tic ;
    % STEP 1. Compute new point
    ObjFun = G(:,y)'*r ;
    ObjFun = ObjFun./Gnorm(y)';
    [maxLOC, indSORT] = max(ObjFun)  ;
    i = y(indSORT(1)) ;
    % STEP 2.  Update alpha and H  (unrestricted least-squares)
    if k==1
        alpha =  G(:,i)\b ;
        H = 1/(G(:,i)'*G(:,i)) ;
    else
        [H alpha] = UpdateWeightsInverse(G(:,z),H,G(:,i),alpha,r) ;
    end
    % STEP 3. Move i from set y to set z
    z = [z;i] ;     y(indSORT(1)) = [] ;
    
    % STEP 4. Find possible negative weights
    n = find(alpha<=0) ;
    if  ~isempty(n)
        % STEP 5
        y = [y; z(n)];  z(n)=[] ;
        H = MultiUpdateInverseHermitian(H,n) ;
        % Recomputing alpha
        alpha = H*(G(:,z)'*b );
        %         if ITER_MINIMUM > NITERATIONS_NO_MORE_POINTS    % before
        %         nov-4th-2021
        %             disp('--------------------------------------------------------------')
        %             disp(['The algorithm cannot proceed with the current set of points'])
        %             disp(['We enlarge the set of candidates with the complementary set'])
        %             y = [y;yCOMPL] ;
        %             ITER_MINIMUM = 0 ;
        %         end
        % ITER_MINIMUM = ITER_MINIMUM + 1;
    else
        % ITER_MINIMUM = 0 ;
    end
    
    if (ITER_MINIMUM > NITERATIONS_NO_MORE_POINTS && COMPLEMENTARY_ALREADY_ADDED ==0 )  || isempty(y)
        disp('--------------------------------------------------------------')
        disp(['The algorithm cannot proceed with the current set of points'])
        disp(['We enlarge the set of candidates with the complementary set'])
        y = [y;yCOMPL] ;
        ITER_MINIMUM = 0 ;
        COMPLEMENTARY_ALREADY_ADDED = 1;
    end
    
    
    
    %     if length(z) > MAX_NUMBER_POINTS  % Removed 4th-Nov-2021
    %         ITER_MINIMUM = 0 ;
    %     else
    %         ITER_MINIMUM = ITER_MINIMUM + 1;
    %         % This happens when a point has been removed (because it was negative)
    %     end
    
    iOLD = i ;
    % STEP 6
    r = b-G(:,z)*alpha ;
    nerror = norm(r)/norm(b) ; % Relative error (using r and b)
    
    if DATA.STORE_INFO_ITERATIONS == 1
        iterLOC = length(z) ;
        DATAOUT.HistoryPoints{iterLOC} =  z  ;
        DATAOUT.HistoryWeights{iterLOC} = alpha.*sqrt(W(z))   ;
        DATAOUT.Error(iterLOC) =  nerror ;
    end
    
    
    % STEP 7
    disp(['k = ',num2str(k),', m=',num2str(length(z)),' ,','; error n(res)/n(b) (%) = ',...
        num2str(nerror*100,10) ]) ;
    ERROR_GLO(k) = nerror ;
    NPOINTS(k) =  length(z) ;
    
    mPOS = length(z) ;
    % MAX_NUMBER_POINTS = max(mPOS,MAX_NUMBER_POINTS) ;
    if  nerror >= nerrorMIN
        ITER_MINIMUM = ITER_MINIMUM + 1 ;
    else
        nerrorMIN = nerror ;
        ITER_MINIMUM = 0 ;
    end
    
    
    
    k = k + 1 ;
    
    %     if length(z) == m
    %         dbstop('88')
    %         disp('')
    %     end
    
end


if  k>= NITERACIONES
    PROPORpoints = length(yORIG)/length(W)*100;
    error(['NO CONVERGENCE. ENLARGE THE SET OF CANDIDATE POINTS  (NOW IT CONTAINS ',num2str(PROPORpoints),' % of the total number of Gauss points)'])
end


w = alpha.*sqrt(W(z)) ;

disp(['Total number of iterations =',num2str(k)])

PLOT_error = 0 ;

if PLOT_error == 1
    
    figure(500)
    hold on
    xlabel('Number of points')
    ylabel('Error (%)')
    plot(NPOINTS,ERROR_GLO*100,'k')
end

DATAOUT.kiteraciones = k ;


% EXAMINE_WEIGHTS_COMPLEMENTARY =1;
% if EXAMINE_WEIGHTS_COMPLEMENTARY ==1
%     [III,JJJ ]= intersect(z,yCOMPL);
%     wCOMPL = w(JJJ) ;
%     meanW = mean(w)
%     meanWcompl = mean(wCOMPL)
% end

% DATA = DefaultField(DATA,'SingularValuesSnapshotMatrix',[]) ;
% if ~isempty(DATA.SingularValuesSnapshotMatrix)
%     % Error pondered by the singular values
%     S = DATA.SingularValuesSnapshotMatrix/DATA.SingularValuesSnapshotMatrix(1) ;
%     if length(S) ~= size(G,1)
%         S = [1; S] ; % The matrix has been expanded with the orthogonal complement
%     end
%     bEX = b.*S ;
%     rEX = r.*S;
%     ErrorEX = norm(rEX)/norm(bEX)*100;
%
%
%
% end

% figure(501)
% hold on
% xlabel('Number of iterations')
% ylabel('Number of points')
% plot(NPOINTS,'k')


