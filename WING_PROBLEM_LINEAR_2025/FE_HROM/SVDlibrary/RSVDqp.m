function [U,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(A,epsilon,DATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a matrix A, RSVDqp returns   a factorization  [U,S,V] of the form
%   A =  U*S*V^T + E   (E = 0 when no truncation is introduced)
% --------------------------
% INPUTS
%------------------------------------------------------------------------------
%  A:    Matrix to be factorized. It can be given as a single numeric array,
%       or as a  p x q cell array containing a conforming partition of A, i.e.:
%
%     A = {A_11, A_12 ... A_1q
%       A_21, A_22 ....A_2q
%       ... .... ..... ...
%       A_p1, A_p2 ... A_pq}
%
%     In turn, each entry A_ij may be a numeric array, or the name of a MAT-file containing
%     the numeric array. The second option is preferred when the whole
%     matrix does not fit into fast memory
%
% -------------------------------------------------------------------------------
%  epsilon :
%     If A is a single numeric array, epsilon is a positive number (below 1) that
%     indicates the error threshold in approximating the whole matrix
%     (in the frobenius norm), i.e.
%     norm(E,'fro')/norm(A,'fro') <= epsilon, where E = A - U*S*V^T
%
%     If A is a cell array, size(epsilon) = size(A), and
%     norm(E{i,j},'fro')/norm(A{i,j},'fro') <= epsilon(i,j)
%
%     NOTES:
%     a) If all the entries of epsilon are zero, then the global error has to
%     be specified by variable  DATA.EPSILON_GLO   (by default = 0)
%
%     b) If DATA.EPSILON_GLO_isRELATIVE =0, the global threshold is given
%     in absolute value
%
% -----------------------------------------------------------------------------
%
%  OUTPUTS
%  ------
%  U --> (M x r)  Matrix of left singular vectors (it approximately spans
%  the column space of A); r denotes the rank of the approximation.
%  V --> (N x r)  Matrix of right singular vectors (it approximately spans
%  the row space of A)
%  S ---> Vector of singular values  (r x 1)
%  ETIME --> Structure containing the times required for each operation
%  eSVD --> Approximation error (Frobenius norm)
%
%  RankMatrix --> rank(Q'*A*P). If epsilon = zeros(size(A)), then
%  RankMatrix = rank(A)
%
%  NOTES: If epsilon(i,j) = 0 , then the obtained
%  factorization is the classical truncated SVD
%  ----------------
%  OPTIONAL INPUTS, (with default values)
%  -----------------
% 1)  DATA.USE_SLOW_MEMORY.STORE_ALSO_SECOND_STEP = 0
% 2)  DATA.Rini = 0
% 3)  DATA.rho_est = 0.05
% 4)  DATA.TRANSPOSE = 0
% 5)  DATA.COMPUTE_V_LAST_SVD = 1
%
%  When each A{i,j} is a .MAT
%  filename,    the contents of A{i,j} are discarded from memory
%  after its corresponding incremental basis vectors are computed.
%  The product of the transponse of basis matrix Qi  times the i-th row
%  block is stored in memory if variable
%  DATA.USE_SLOW_MEMORY.STORE_ALSO_SECOND_STEP = 1 ;
%
%    DATA.Rini = Initial estimation for the rank of submatrix A{1,1}. The
%    default value is min(DATA.rho_est*ceil(max(size(A{1,1})))).
%
%  DATA.TRANSPOSE -> The factorization is calculated on the transpose of
%  the matrix
%
%  DATA.COMPUTE_V_LAST_SVD -->  If 0, the matrix of right singular vectors
%  is not calculated
%
%
%
%  Written by Joaquín A. Hernández Ortega, Dec. 2016
%  UPC/CIMNE. jhortega@cimne.upc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    load('tmp1.mat')
    
end

if exist('DefaultField') == 0
    addpath('AUXFUN')
end


if nargin == 2
    DATA = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT INPUTS
% ----------------------------------------------------
DATA = DefaultField(DATA,'EPSILON_GLO',0) ; % See preabmble
DATA = DefaultField(DATA,'USE_SLOW_MEMORY',[]) ;  %
DATA.USE_SLOW_MEMORY = DefaultField(DATA.USE_SLOW_MEMORY,'STORE_ALSO_SECOND_STEP',0) ;  %
DATA = DefaultField(DATA,'EPSILON_GLO_isRELATIVE',1) ; % For multi-level version
DATA = DefaultField(DATA,'COMPUTE_V_LAST_SVD',1) ; % It computes the matrix of right-singular vectors
DATA = DefaultField(DATA,'TWOLEVELMETHOD',[]) ; % Two-level method
DATA.TWOLEVELMETHOD = DefaultField(DATA.TWOLEVELMETHOD,'ACTIVE',0) ; % Two-level method
DATA = DefaultField(DATA,'Rini',0) ; % Initial estimation for the rank
DATA = DefaultField(DATA,'TRANSPOSE',0) ; %
DATA = DefaultField(DATA,'PARALLEL',0) ; %
DATA = DefaultField(DATA,'rho_est',0.05) ; %
DATA = DefaultField(DATA,'WHOLE_DECOMPOSITION_GLO',0) ; %
DATA = DefaultField(DATA,'alphaA',[]) ; %
DATA = DefaultField(DATA,'betaA',[]) ; %
DATA = DefaultField(DATA,'SetGammaToZero',0) ; %
DATA = DefaultField(DATA,'normAoriginal',[]) ; %
DATA = DefaultField(DATA,'HIDE_OUTPUT',1) ; %




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETIME.TOTAL = tic ;
% Preprocessing input data (check whether it is a cell array, or a numeric array and other operations)
% -------------------------------------------------------------------------------------------------------
NAMEF = [] ;
DATA.USE_SLOW_MEMORY.ACTIVE = 0 ;
if iscell(A)
    CHECK = cellfun(@ischar,A) ;
    if  any(CHECK)
        [ aaa bbb]= find(CHECK==1) ;
        NAMEF = A{aaa(1),bbb(1)} ;
        DATA.USE_SLOW_MEMORY.ACTIVE = 1 ;
    end
end
if isnumeric(A) || (iscell(A) && isempty(NAMEF) && size(A,1)==1 & size(A,2)==1 )
    if  size(A,1)==1 & size(A,2)==1
        A = A{1,1} ;
    end
    DATA.PATH_STORE_AUX  = [] ;
    p = 1; q = 1;
    DATA.TRANSPOSE = 0 ;
else
    [p,q] = size(A) ;
    if ~isempty(NAMEF)
        NAMEF = fileparts(NAMEF) ;
        DATA.PATH_STORE_A = NAMEF ;
        DATA.PATH_STORE_AUX = [DATA.PATH_STORE_A,'AUXtmpj/'] ;  % Path in which temporary files are stored
    else
        DATA.PATH_STORE_AUX  = [] ;
        
    end
    % dbstop('34')
    if DATA.TRANSPOSE == 1
        p_old = p ; q_old = q ;
        p = q_old ; q = p_old ;
        A  = A' ;
        epsilon = epsilon' ;
    end
end
% Creating folder for storing temporary files
if ~isempty( DATA.PATH_STORE_AUX ) && ~exist(DATA.PATH_STORE_AUX)
    mkdir(DATA.PATH_STORE_AUX)
end


[alpha beta] = DetermineAlphaBeta(A) ;
DATA.alphaA = alpha ;
DATA.betaA = beta ;



%-----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = [] ; P = [] ; alpha_c = [] ; ETIME_Qoper =[] ; ETIME_Poper = [] ;
if p >1 | q >1
    % while iter <= DATA.ITERMAX & SIZEA>=DATA.MaxSizeMatrixGb
    if DATA.HIDE_OUTPUT == 0
        disp('----------------')
        disp('COMPUTATION OF Q and Q^T*A  (basis for column space of A)')
        disp('----------------')
    end
    tloc = tic ;
    [Q, A, gamma, alpha_1, beta_1,  Rmax, ETIME_Qoper] ...
        = RORTHmatrix(A,epsilon,DATA) ;
    
    ETIME.Q = toc( tloc) ;
    if DATA.HIDE_OUTPUT == 0
        disp('----------------')
        disp('COMPUTATION OF P and Q^T*A*P (basis for the row space of A)')
        disp('----------------')
    end
    DATA.USE_SLOW_MEMORY.STORE_ALSO_SECOND_STEP = 0;
    DATA.EPSILON_ABSOLUTE = 1 ;
    tloc = tic ;
    DATA.Rini = Rmax ;
    % dbstop('148')
    %DATA.TryToPrecomputeQtB = 0;
    if DATA.SetGammaToZero == 1
        gamma = zeros(size(gamma)) ;
    end
    
    TRANSPOSE = DATA.TRANSPOSE ; DATA.TRANSPOSE = 0 ;  %DATA.PARALLEL = 0 ;
    DATA.alpha = alpha_1 ; DATA.beta = beta_1 ;
    [P A gamma_c alpha_c beta_c Rmax ETIME_Poper]  =  RORTHmatrix(A,gamma',DATA) ;
    DATA.TRANSPOSE  = TRANSPOSE;
    ETIME.P = toc( tloc) ;
    
    [M,~] = cellfun(@size,Q) ;  M  = sum(M) ;
    [N,~] = cellfun(@size,P) ;  N = sum(N) ;
    
else
    if ~isnumeric(A) & ~isnumeric(A{1,1})
        SSS = load(A{1,1}) ;
        fff = fieldnames(SSS) ;
        A = SSS.(fff{1}) ;
    end
    DATA.EPSILON_GLO = epsilon(1); % Modification 16-Jan-2021 !!!!!
    [M  N]= size(A) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DATA.HIDE_OUTPUT == 0
    disp('----------------')
    disp('SVD of A <-- Q^T A P')
    disp('----------------')
end
eLOC = DATA.EPSILON_GLO ;
if ~isnumeric(A)
    A = cell2mat(A) ;
end
tloc = tic ;
c = norm(A,'fro') ;
%dbstop('159')
mu = (max([M N])*eps(c))  ;  % Machine epsilon parameter
if DATA.EPSILON_GLO_isRELATIVE == 1
    if  eLOC == 0
        e0 = mu ;
    else
        e0 = c*eLOC ;
    end
else
    e0 = eLOC ;
end
% -----------------------------------------
% Estimation rank A
if DATA.Rini == 0
    Rest = ceil(DATA.rho_est*min(size(A))) ;
end
if ~isempty(alpha_c)
    Rest = max([alpha_c beta_c]) ;
else
    Rest = min(Rest,max(size(A))) ;
end

DATA.COMPUTE_V_SVD = DATA.COMPUTE_V_LAST_SVD;


DATA.WHOLE_DECOMPOSITION = DATA.WHOLE_DECOMPOSITION_GLO ;

%[U,S,V,eSVD,RankMatrix] = RSVDt(A,e0,mu,Rest,DATA) ;
[U,S,V,eSVD,RankMatrix] = RSVDT(A,e0,mu,Rest,DATA) ;



% end
if DATA.HIDE_OUTPUT == 0
    
    disp('---------------------------------------')
    disp('DONE ... ')
    disp('--------------------------------------')
end
ETIME.QAP = toc(tloc);
if   DATA.HIDE_OUTPUT == 0
    
    disp('Computing left singular vectors U = X*Ubar')
end
tloc = tic ;
if ~isempty(Q)
    if iscell(Q)
        U = QUprod(Q,U,DATA.alphaA) ;
    else
        U = Q*U    ;
    end
end

if  DATA.HIDE_OUTPUT == 0
    
    disp('Done ')
end
ETIME.U = toc(tloc) ;
tloc = tic ;
if  DATA.COMPUTE_V_LAST_SVD == 1
    if   DATA.HIDE_OUTPUT == 0
        
        disp('---------------------------------------')
        disp('Computing right singular vectors V = P*Vbar')
    end
    if ~isempty(P)
        if iscell(P)
            V = QUprod(P,V,DATA.betaA) ;
        else
            V = P*V    ;
        end
    end
    disp('Done ')
else
    P = [] ;
end
ETIME.V = toc(tloc) ;

if  DATA.TWOLEVELMETHOD.ACTIVE == 0 & ~isempty(DATA.PATH_STORE_AUX)
    rmdir(DATA.PATH_STORE_AUX,'s')
end

ETIME.TOTAL = toc(ETIME.TOTAL) ;
ETIME.Qoper= ETIME_Qoper ;
ETIME.Poper= ETIME_Poper ;

if DATA.TRANSPOSE == 1
    Unew = V ;
    V = U ;
    U = Unew ;
end

DATAOUT.mu = mu ;




