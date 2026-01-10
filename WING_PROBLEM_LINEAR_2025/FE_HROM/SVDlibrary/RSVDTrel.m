function  [U,S,V,eSVD,Rsup] = RSVDTrel(A,e0,DATA)
% RSVDT computes the truncated singular value decomposition (SVD) of matrix
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

 
% DEFAULT INPUT DATA
%---------------------
if nargin == 0
    load('tmp1.mat') 
elseif nargin == 1
    e0 = 0 ; 
    DATA = [] ; 
elseif nargin == 2
    DATA = [] ; 
%         DATA = [] ;
% elseif nargin == 1
%     e0 = 0 ;     DATA = [] ;
% elseif nargin == 2
%     mu = 0 ;  R = 0 ;    DATA = [] ;
% elseif nargin == 3
%     R =0 ;   DATA= [] ;
% elseif nargin == 4
%     DATA= [] ;
end

 mu=0 ; R = 0 ;

DATA = DefaultField(DATA,'COMPUTE_V_SVD',1) ;
DATA = DefaultField(DATA,'RELATIVE_SVD',1) ;
DATA = DefaultField(DATA,'USE_ALWAYS_RANDOMIZATION',1) ; 
% END DEFAULT INPUT DATA
% ------------------------
Rsup = min(size(A)) ;
if R >= Rsup  && DATA.USE_ALWAYS_RANDOMIZATION == 0
    Q = 'full' ; %    
else
    if R>=Rsup
        % This means that the algorithm should take all possible modes. 
        R = 0 ;
    end    
    [Q B eORTH a]= RORTH(A,mu,R,DATA) ; % Randomized orthogonalization (machine precision parameter = mu)
    Rsup = size(Q,2) ;
end

if isempty(Q)
    U = [] ;     S = [] ;     V = [] ;
    eSVD = 0 ;
else
    if  isnumeric(Q)
        if DATA.RELATIVE_SVD ==1 && e0 >0 %
            e0 = e0*a ;
        end
        if e0 > eORTH;   e = sqrt(e0^2 - eORTH^2) ;   else   e = e0 ; end
        
        %    dbstop('95')
        DATA.RELATIVE_SVD = 0 ;
        DATA.SVD.MaxSize = max(size(A))  ;
        [U,S,V,eSVDb] = SVDT(B,e,DATA);
        U = Q*U ;
        eSVD = sqrt(eORTH^2 + eSVDb^2) ;
    else
        % A appears to be full rank
        [U,S,V,eSVD] = SVDT(A,e0,DATA);
        
    end
end


