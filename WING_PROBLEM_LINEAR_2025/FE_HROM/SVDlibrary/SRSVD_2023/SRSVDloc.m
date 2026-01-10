function [U,S,V] = SRSVDloc(A,RELTOL,DATA)
%--------------------------------------------------------------------------
% function [U, S, V] = SRSVDloc(A, RELTOL, DATA)
%
% PURPOSE:
%   Computes a Structured Randomized Singular Value Decomposition (SRSVD)
%   of a matrix **partitioned by columns**, i.e. A = {A₁, A₂, ..., A_q}, 
%   where each Aᵢ is a submatrix with the same number of rows.
%
%   This is a **memory-aware** algorithm designed for very large matrices
%   (e.g., those that do not fit entirely in memory), and is part of the
%   block-wise randomized SVD framework described in Appendix A.3 of:
%
%      Hernández et al., 2024, *Continuous Empirical Cubature Method (CECM)*
%
%   The procedure:
%     1. Computes an orthonormal basis Q of the column space of A using
%        a streaming Gram-Schmidt variant (see `SRORTH`).
%     2. Projects A onto the Q subspace to form a compressed matrix L.
%     3. Applies truncated SVD to L with tolerance derived from RELTOL.
%     4. Reconstructs the dominant left singular vectors: U = Q * Ũ.
%
% INPUTS:
%   - A : Cell array {A₁, ..., A_q}, where each Aᵢ is a numeric matrix
%         with size (n_rows × n_colsᵢ). The union of Aᵢs spans the full
%         matrix A to be factorized.
%
%   - RELTOL : Relative truncation threshold
%         · If RELTOL is a scalar, global Frobenius-norm-based truncation is applied:
%               ‖A - U*S*Vᵗ‖_F ≤ RELTOL * ‖A‖_F
%         · If RELTOL is a vector of length q, blockwise local truncation
%           is applied in the orthonormalization (used in `SRORTH`).
%
%   - DATA : Struct with optional settings, including:
%         · DATA.HIDE_OUTPUT : (default = 0) Suppresses console output.
%         · DATA.SVD.MaxSize : (optional) Matrix size for adaptive tolerance.
%
% OUTPUTS:
%   - U : (n_rows × r) Left singular vectors (approximate column space of A)
%   - S : (r × 1) Vector of dominant singular values
%   - V : (n_cols × r) Right singular vectors
%
% METHOD:
%   1. **Orthonormalization via SRORTH**:
%         Computes Q such that Qᵗ * A ≈ L
%         L is a compact matrix (cell of projected blocks)
%         γ measures the orthogonalization residuals
%
%   2. **SVD on the compressed matrix**:
%         Performs truncated SVD on cell2mat(L) using `SVDTloc`,
%         with adjusted tolerance based on Frobenius norm and orthogonalization error.
%
%   3. **Reconstruction**:
%         Left singular vectors of A are recovered as Q * Ũ,
%         where Ũ are the left singular vectors of L.
%
% REMARK:
%   This function is optimized for high-dimensional data and reduced-order modeling
%   in nonlinear problems (e.g., hyperreduction), where A represents a basis of
%   candidate modes or snapshots to be filtered via randomized projections.
%
% SEE ALSO:
%   - SRORTH      (streaming orthonormalization)
%   - SVDTloc     (truncated SVD for block-matrix input)
%   - SRSVD       (general wrapper for structured randomized SVD)
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, Oct. 2022 – UPC/CIMNE – jhortega@cimne.upc.edu
%--------------------------------------------------------------------------

% See comments in RSVDTrowblock.m 

if nargin == 0
    load('tmp.mat')
end
if nargin == 2
    DATA = [] ;
end

if length(RELTOL) == length(A)
    % Block-wise tolerance
    TOL_BLOCK = RELTOL ;
    TOL_GLO = 0 ;
else
    TOL_BLOCK = 0*ones(size(A)) ;
    TOL_GLO = RELTOL ;
end

DATA = DefaultField(DATA,'HIDE_OUTPUT',0) ; % = 0;
[Q,L,gamma,NormA,ETIME] = SRORTH(A,TOL_BLOCK,DATA) ;
% Q is an "exact" orthogonal basis matrix for the column space of A  y
% TOL_BLOCK =[0,0 ...0]; otherwise, it is an approximated basis matrix,
% with approximation given by TOL_BLOCK
NormA = norm(NormA) ;


e0 = TOL_GLO ;
eORTH = norm(gamma) ;
%if  size(Q,2) < ncolsTOT
%  if DATA2.RELATIVE_SVD ==1 && TOL_GLO >0 %
e0 = TOL_GLO*NormA;
% end
if e0 > eORTH;   e = sqrt(e0^2 - eORTH^2) ;   else   e = e0 ; end

%    dbstop('95')
DATA.RELATIVE_SVD = 0 ;
[sA,sB] = cellfun(@size,A,'UniformOutput',false) ;
nrowsTOT = sA{1} ;
ncolsTOT = sum(cell2mat(sB)) ;


DATA.SVD.MaxSize = max([nrowsTOT,ncolsTOT])  ;
[U,S,V,eSVDb] = SVDTloc(cell2mat(L),e,DATA);
U = Q*U ;
eSVD = sqrt(eORTH^2 + eSVDb^2) ;
%else
% A appears to be full rank
%   [U,S,V,eSVD] = SVDT(A,e0,DATA);

%end
