function [U,S,V,eSVD] = SVDTloc(B,epsilon,DATA)
%--------------------------------------------------------------------------
% function [U, S, V, eSVD] = SVDTloc(B, epsilon, DATA)
%
% PURPOSE:
%   Computes a **truncated Singular Value Decomposition** (SVD) of matrix B,
%   retaining only the dominant singular values and associated singular vectors
%   according to an absolute or relative tolerance threshold.
%
%   This function is part of the CECM framework and is used in hyperreduction
%   contexts to compress a dense basis or snapshot matrix into a low-rank
%   representation with controllable error.
%
% INPUTS:
%   - B : Matrix (M×N) to be factorized.
%   - epsilon : Tolerance controlling truncation of the SVD.
%         If epsilon = 0 → uses machine precision to determine rank.
%         If epsilon > 0 → truncates based on Frobenius norm.
%
%   - DATA : Struct of options with fields:
%         · COMPUTE_U : 1 (default) → compute left singular vectors U
%         · COMPUTE_V : 1 (default) → compute right singular vectors V
%         · RELATIVE_SVD : 0/1 → if 1, tolerance is relative to ‖B‖_F
%         · SVD.MaxSize : scalar, used for defining rank tolerance
%
% OUTPUTS:
%   - U : (M × K) matrix with K left singular vectors (if COMPUTE_U = 1)
%   - S : (K × 1) vector of dominant singular values
%   - V : (N × K) matrix with K right singular vectors (if COMPUTE_V = 1)
%   - eSVD : Truncation error estimate (Frobenius norm of discarded singular values)
%
% METHOD OVERVIEW:
%   1. Performs thin SVD on B or Bᵗ depending on its shape (tall vs wide).
%   2. Converts singular values matrix to vector form.
%   3. Computes numerical rank R based on machine epsilon.
%   4. Applies truncation based on epsilon (absolute or relative).
%   5. Truncates U, S, and V to reduced-rank form.
%
% DETAILS:
%   - The truncation condition is derived from:
%       √(∑_{i=K+1}^N σᵢ²) ≤ epsilon
%     where σᵢ are singular values.
%
%   - If RELATIVE_SVD = 1, epsilon is scaled by ‖B‖_F.
%
% REFERENCES:
%   - Hernández et al., 2024, *Continuous Empirical Cubature Method (CECM)*.
%     See Appendix A.2 and A.3 for theory on rank-adaptive compression.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, March 2017
%   UPC/CIMNE, Barcelona — jhortega@cimne.upc.edu
%--------------------------------------------------------------------------
% Comments by ChatGPT4


%  SVDT(B,epsilon) computes a truncated SVD of B (with   tolerance epsilon)
%  It employs the built-in matlab function "svd", but it only returns the
% singular values  truncated according to the   the specified tolerance
%
% S = SVDT(B,epsilon) gives the truncated vector of singular values
%
% [U,S] = SVDT(B,epsilon) furnishes the truncated matrix of left singular vectors and
% the associated singular values
%
% [U,S,V] = SVDT(B,epsilon)  returns the full factorization (truncated).
%
% By default,  norm(U*diag(S)*V'-B,'fro') <= epsilon
% If epsilon = 0, then the truncation criterion is the epsilon machine
% parameter
% If the third argument REL = 1, (by default REL = 0), then
% norm(U*diag(S)*V'-B,'fro') <= epsilon*norm(B,'fro')*
%
% -----------------------------------
% Written by Joaquín  A. Hernández, March 2017.
% UPC/CIMNE, Barcelona, Spain
% jhortega@cimne.upc.edu
if nargin == 0
    load('tmp1.mat')
    epsilon = 0 ;  DATA = [] ;
end
if nargin ==1
    epsilon = 0 ;  DATA = [] ;
elseif nargin ==2
    DATA = [] ;
end

DATA = DefaultField(DATA,'COMPUTE_U',1) ;
DATA = DefaultField(DATA,'COMPUTE_V',1) ;
DATA = DefaultField(DATA,'RELATIVE_SVD',0) ;
DATA = DefaultField(DATA,'SVD',[]) ;
DATA.SVD = DefaultField(DATA.SVD,'MaxSize',max(size(B))) ;



CALCU =DATA.COMPUTE_U ; CALCV = DATA.COMPUTE_V ; U = [] ; V= [] ;
dimMATRIX = DATA.SVD.MaxSize;
M = size(B,1); N = size(B,2);
if M>=N
    if CALCU==1 & CALCV==1
        [U, S,V] = svd(B,'econ') ;  % U --> M xN, V --> N x N
    elseif CALCU==1 & CALCV==0
        [U, S] = svd(B,'econ') ;  % U --> M xN, V --> N x N
    else
        S = svd(B,'econ') ;
    end
else
    % If N>M, it proves more efficient to perform the SVD of B'
    if CALCU==1 & CALCV==1
        [V, S,U] = svd(B','econ') ; % U --> M x M, V --> Nx M
    elseif  CALCU==1 & CALCV==0
        [~, S,U] = svd(B','econ') ; % U --> M x M, V --> Nx M
    else
        S = svd(B','econ') ;
    end
end
if size(S,2)>1
    S = diag(S) ; % S = [S1 S2 ... SR]'
end

tol = dimMATRIX*eps(max(S));

R = sum(S > tol);  % Definition of numerical rank

eSVD = 0 ;

%dbstop('73')
if  epsilon ==0
    K = R ;
else
    %    dbstop('70')
    SingVsq =  (S.*S) ;
    SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
    normEf2 = sqrt(cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
    
    if DATA.RELATIVE_SVD ==1
        epsilon = epsilon*normEf2(end) ;
        
    end
    T = (sum(normEf2<epsilon)) ;
    K = length(S)-T ;
    
    
end
K = min(R,K);

%dbstop('88')
if length(S) >K
    eSVD = sqrt(sum(S(K+1:end).^2)) ;
end

S = S(1:K);
if  ~isempty(U) & isempty(V)
    U = U(:,1:K)  ;
    %  varargout={U,S};
elseif ~isempty(U) & ~isempty(V)
    U = U(:,1:K)  ;
    V =  V(:,1:K) ;
    % %
    % else
    %     varargout = {S} ;
end




