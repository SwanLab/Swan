function [U,S,V ,Mchol]= WSVDT(A,M,DATALOC)
%--------------------------------------------------------------------------
% function [U, S, V, Mchol] = WSVDT(A, M, DATALOC)
%
% PURPOSE:
%   Computes the **Weighted Singular Value Decomposition (WSVD)** of matrix A,
%   using a symmetric positive-definite weight matrix M to induce a new inner product.
%
%   That is, this decomposition ensures the columns of U are **M-orthonormal**:
%     Uᵗ M U = I
%
% USAGE:
%   - [U, S, V] = WSVDT(A, M)
%       Computes the M-weighted SVD of A.
%
%   - [U, S, V, Mchol] = WSVDT(A, M, DATALOC)
%       Also returns the Cholesky factor of M used in the transformation.
%
% INPUTS:
%   - A       : [m x n] Input matrix to decompose.
%   - M       : [m x m] Symmetric positive-definite matrix defining the weight.
%   - DATALOC : (Optional) Struct with additional options:
%       * DATALOC.TOL           : Truncation tolerance for SVDT (default = 0)
%       * DATALOC.RELATIVE_SVD  : Use relative tolerance (default = 1)
%       * DATALOC.Mchol         : If provided, uses this Cholesky factor of M
%
% OUTPUTS:
%   - U      : Matrix with M-orthonormal left singular vectors (Uᵗ M U = I)
%   - S      : Vector of singular values (possibly truncated)
%   - V      : Matrix with right singular vectors (Vᵗ V = I)
%   - Mchol  : Cholesky factor of M (such that M = Mcholᵗ * Mchol)
%
% REMARKS:
%   - If M is empty, standard Euclidean SVD is performed.
%   - Internally transforms A → Ā = L * A with M = Lᵗ * L, then applies standard SVD.
%   - The transformation is undone at the end: U = L⁻¹ * Ū.
%
% SEE ALSO:
%   - SVDT.m (used internally for truncated SVD)
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE (Barcelona, Spain)
%   April 2024
%   jhortega@cimne.upc.edu
%--------------------------------------------------------------------------

if nargin == 2
    DATALOC = [] ;
end
% Weighted SVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mutually  M-orthogonal
NO_CHOL = 0 ; 
DATALOC = DefaultField(DATALOC,'Mchol',[]) ;
if  isempty(DATALOC.Mchol)
    if ~isempty(M)
    disp(['Computing Choleski Decomposition ....'])
    Mchol = chol(M) ;
    else
        disp(['  Choleski Decomposition  not computed... (= Euclidean) ....'])
        Mchol = 1; 
        NO_CHOL = 1; 
    end
else
    Mchol = DATALOC.Mchol ;
    if isempty(Mchol)
         disp(['  Choleski Decomposition  not computed... (= Euclidean) ....'])
        Mchol = 1; 
        NO_CHOL = 1 ; 
    end 
end
%Abar = Mchol*A ;
DATALOC = DefaultField(DATALOC,'TOL',0) ;
TOL = DATALOC.TOL ;
DATALOC = DefaultField(DATALOC,'RELATIVE_SVD',1) ;

[U,S,V] = SVDT( Mchol*A,TOL,DATALOC) ;
if  NO_CHOL == 0  % JAHO, 22-Apr-2204
U = Mchol\U ;
end

% Sjump  = zeros(size(S));
% Sjump(1) = 1;
% Sjump(2:end) = S(2:end)./S(1:end-1) ;
