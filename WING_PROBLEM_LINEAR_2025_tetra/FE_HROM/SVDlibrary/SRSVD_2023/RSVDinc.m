function  [U,S,V,eSVD,Rsup] = RSVDinc(A,e0,mu,R,DATA)
%--------------------------------------------------------------------------
% function [U,S,V,eSVD,Rsup] = RSVDinc(A,e0,mu,R,DATA)
%
% PURPOSE:
%   Computes a **low-rank approximation** of matrix A using an *incremental randomized SVD*:
%
%       A ≈ U * diag(S) * Vᵗ      such that ‖A - USVᵗ‖_F ≤ e0
%
%   This is a key subroutine in the *CECM* methodology for obtaining fast and
%   memory-efficient spectral decompositions with guaranteed accuracy.
%
%   The method first builds an orthonormal basis `Q` for the column space of A
%   via `RORTHinc`, and then computes the SVD of the projected matrix `B = QᵗA`.
%   If truncation threshold `e0` is provided, only the minimal number of modes
%   needed to satisfy ‖A - USVᵗ‖ ≤ e0 are retained.
%
% INPUT:
%   - A : [M × N] matrix to approximate
%
%   - e0 : Desired upper bound for ‖A - USVᵗ‖_F
%         If 0, it will be set to mu (machine epsilon)
%
%   - mu : Machine precision estimator, used to set default e0 if needed.
%         If empty, computed as: mu = min(size(A)) * eps(norm(A,'fro'))
%
%   - R : Estimated upper bound on the numerical rank of A.
%         If 0 or missing, initialized as ceil(ρ * min(size(A)))
%
%   - DATA : Structure with fields:
%         * DATA.rho                : (default 0.05) used when R = 0
%         * DATA.NITER              : max # iterations for orthonormalization
%         * DATA.COMPUTE_V_SVD      : flag to compute right singular vectors
%         * DATA.RELATIVE_SVD       : if 1, e0 is interpreted as relative tolerance
%         * DATA.USE_ALWAYS_RANDOMIZATION : if 0, disable randomized scheme when rank is full
%
% OUTPUT:
%   - U : [M × r] matrix of truncated left singular vectors
%   - S : [r × 1] vector of retained singular values
%   - V : [N × r] matrix of truncated right singular vectors (if requested)
%   - eSVD : Frobenius norm of the residual error ‖A - USVᵗ‖_F
%   - Rsup : Effective rank used in the randomized projection (cols of Q)
%
% METHOD OVERVIEW:
%   1. Compute orthonormal basis Q via randomized orthogonalization (RORTHinc)
%   2. Project A to B = QᵗA
%   3. Compute low-rank SVD of B, truncate according to e0
%   4. Map back U = Q * U_B
%
% REMARKS:
%   - If A is numerically full-rank and DATA.USE_ALWAYS_RANDOMIZATION = 0, a
%     full (non-randomized) truncated SVD is used instead.
%
%   - When DATA.RELATIVE_SVD = 1, the user-provided e0 is interpreted as a 
%     relative tolerance with respect to ‖A‖_F
%
%   - Used internally in: SRORTH, SRSVDloc, and other CECM routines
%
% REFERENCES:
%   - Hernández et al., 2024. "CECM: A Continuous Empirical Cubature Method..."
%     Appendix A.2 (Randomized SVD with controlled truncation)
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, Dec 2016 / Rev. 2022
%   UPC / CIMNE, Barcelona
%   jhortega@cimne.upc.edu
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------------------------


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
% Written by Joaquín  A. Hernández, December 2016/2022.
% UPC/CIMNE, Barcelona, Spain
% jhortega@cimne.upc.edu
% ------------------------------------------------------------------

 
% DEFAULT INPUT DATA
%---------------------
if nargin == 0
    load('tmp1.mat') 
        DATA = [] ;
elseif nargin == 1
    e0 = 0 ;  mu=0 ; R = 0 ;    DATA = [] ;
elseif nargin == 2
    mu = 0 ;  R = 0 ;    DATA = [] ;
elseif nargin == 3
    R =0 ;   DATA= [] ;
elseif nargin == 4
    DATA= [] ;
end

DATA = DefaultField(DATA,'HIDE_OUTPUT',1) ;
DATA = DefaultField(DATA,'COMPUTE_V_SVD',1) ;
DATA = DefaultField(DATA,'RELATIVE_SVD',0) ;
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
    [Q, B, eORTH, a]= RORTHinc(A,mu,R,DATA) ; % Randomized orthogonalization (machine precision parameter = mu)
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
        [U,S,V,eSVDb] = SVDTloc(B,e,DATA);
        U = Q*U ;
        eSVD = sqrt(eORTH^2 + eSVDb^2) ;
    else
        % A appears to be full rank
        [U,S,V,eSVD] = SVDTloc(A,e0,DATA);
        
    end
end


