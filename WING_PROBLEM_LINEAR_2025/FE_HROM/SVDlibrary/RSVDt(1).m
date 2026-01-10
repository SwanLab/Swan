function  [U,S,V,e_svd,RankMatrix] = RSVDt(C,e0,mu,R,DATA)
% RSVDt computes the truncated singular value decomposition (SVD) of matrix
% %
%
% C = U*diag(S)*V' + E(e0,mu)
%
% INPUTS (mandatory)
% ------
% C --> M x N matrix whose factorization is to be calculated
%
%  Optional inputs (with default values)
%  ---------------
%
%  1) e0 = 0
%  2) mu = []
%  3) R  = 0
%  4) DATA.NITER = 10
%  5) DATA.rho = 0.05
%  6) DATA.COMPUTE_V_SVD = 1
%
%  e0 --> Truncation threshold  (only the singular values and associated vectors that renders
%  the frobenious norm of the residual E below this threshold are included in the factorization)
%
%     If e0 = 0, we make e0 = mu (mu is defined below)
%
%  mu --> Machine precision parameter
%        If mu = [] , mu =  min(size(C))*eps(norm(C,'fro'))
%
%  R is an estimation for an upper bound of rank(C) %
%      if nargin == 1 or R=0,  then R = ceil(DATA.rho*min(size(C))),
%
%  DATA.NITER = Maximum number of iterations (rank revealing algorithm).
%
%  DATA.COMPUTE_V_SVD --> Compute matrix of right singular vectors
%
% OUTPUT:
% ------
% U --> (truncated)  matrix of left singular vectors
% S ---> (truncated) vector of singular values
% V --> (truncated) matrix of right singular vectors
% e_svd --> Approximation error (frobenius norm)
% RankMatrix --> Rank matrix.
% ----------------------------------------------------------------
% Written by Joaquín  A. Hernández, December 2016.
% UPC/CIMNE, Barcelona, Spain
% jhortega@cimne.upc.edu
% ------------------------------------------------------------------

%dbstop('50')
if nargin ==0
    load('tmp1.mat')
end

tic
if nargin == 1
    e0 = 0 ;
    mu=[] ;
    R = 0 ;
    DATA = [] ;
elseif nargin == 2
    mu = [] ;
    R = 0 ;
    DATA = [] ;
elseif nargin == 3
    R = 0 ; 
    DATA= [] ;
    elseif nargin == 4
        DATA = [] ; 
end

DATA = DefaultField(DATA,'COMPUTE_V_SVD',1) ;
DATA = DefaultField(DATA,'NITER_RORTH',10) ;
DATA = DefaultField(DATA,'WHOLE_DECOMPOSITION',0) ;


if R >= min(size(C))
    Q = 'full' ; % This has to be revised (unequal submatrices)
else
    [Q D]= RANDQB_adapt(C,R,mu,DATA) ;
end

if isempty(Q)
    U = [] ;     S = [] ;     V = [] ;
    e_svd = 0 ;  RankMatrix = 0 ;
else
    if  isnumeric(Q)
        if isempty(D)
            D = Q'*C ;
        end
        if DATA.COMPUTE_V_SVD == 0
            [U,S] = SVD(D,0);
            V = [] ;
        else
            [U,S,V] = SVD(D,0); %
        end
        U = Q*U ;
    else
        % C appears to be full rank
        if DATA.COMPUTE_V_SVD == 0
            [U,S] = SVD(C,0); %
            V = [] ;
        else
            [U,S,V] = SVD(C,0); %
        end
    end
    %% Truncation
    RankMatrix = length(S);
    if  e0<=mu
        R = length(find(S>mu)) ;
        U = U(:,1:R) ;
        S = S(1:R) ;
        if ~isempty(V)
        V = V(:,1:R) ;
        end
    else
        disp(['Rank matrix =',num2str(RankMatrix)]) ;
        SingVsq =  (S.*S) ;
        SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
        normEf2 = sqrt(cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
        tol = e0 ;
        if tol<=mu
            R = length(S) ;
        else
            T = (sum(normEf2<tol)) ;
            R = length(S)-T ;
        end
    end
    % Actual error
    e_svd = sqrt(sum(S(R+1:end).^2)) ;
    
    if DATA.WHOLE_DECOMPOSITION==0
        U = U(:,1:R);
        S = S(1:R) ;
        if ~isempty(V)
            V = V(:,1:R) ;
        end
    else
        
        RankMatrix = R ;
        
    end
    disp(['Truncated rank=',num2str(R)]) ;
end


