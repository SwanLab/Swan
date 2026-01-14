function [U,S,V] = SRSVD(A,RELTOL,DATA)
%--------------------------------------------------------------------------
% function [U, S, V] = SRSVD(A, RELTOL, DATA)
%
% PURPOSE:
%   Computes a **Sequential Randomized Singular Value Decomposition (SRSVD)** 
%   of a matrix A (possibly partitioned into blocks), such that:
%
%       A ≈ U * diag(S) * V^T     with ||A - U*S*V^T||_F / ||A||_F ≤ RELTOL
%
%   The algorithm is suited for large-scale applications and memory-efficient
%   implementations. It supports both in-memory arrays and disk-stored matrix blocks
%   (e.g., via MAT-files), as described in Appendix A.3 of the referenced paper.
%
% INPUTS:
%   - A : Matrix or 1×q cell array {A₁, A₂, ..., A_q}
%         representing a conforming partition of the matrix A. Each Aᵢ may be:
%         · a numeric matrix (in RAM), or
%         · a string pointing to a MAT-file containing the numeric matrix Aᵢ.
%
%   - RELTOL : Relative truncation threshold (scalar or vector)
%         If scalar: global accuracy control for the Frobenius norm.
%             ||A - U*S*V^T||_F / ||A||_F ≤ RELTOL
%         If vector: block-wise accuracy control: for each block Aᵢ,
%             ||Aᵢ - Uᵢ*Sᵢ*Vᵢ^T||_F / ||Aᵢ||_F ≤ RELTOL(i)
%
%   - DATA : Struct of optional parameters. Relevant fields include:
%         · DATA.RELATIVE_SVD : (default = 1) Controls relative truncation.
%         · DATA.WeightsPremultipy_matrix : optional weights (e.g., √W) to pre-multiply A.
%
% OUTPUTS:
%   - U : Matrix with orthonormal columns spanning the dominant left singular
%         subspace of A (dimension M × r, where r is the numerical rank).
%
%   - S : Vector of singular values corresponding to the retained rank.
%
%   - V : Matrix with orthonormal columns spanning the dominant right singular
%         subspace of A (dimension N × r).
%
% ALGORITHM:
%   The procedure selects between:
%     · Direct randomized SVD (`RSVDinc`) for unpartitioned input.
%     · Blockwise randomized SVD (`SRSVDloc`) when A is given as a block partition.
%
%   In both cases, the decomposition relies on randomized projections and
%   iterative enrichment (see Appendix A.1 and A.3), followed by:
%     - orthonormalization,
%     - computation of the projected SVD,
%     - and optional truncation based on Frobenius-norm-based error control.
%
% REMARKS:
%   · If A is too large for memory, the partitioned format {A₁, ..., A_q}
%     allows streaming evaluation and memory-efficient processing.
%
%   · This function is part of the reduced-order modeling framework introduced in:
%     Hernández et al., *Continuous Empirical Cubature for Nonlinear Model Reduction*
%     (CECM), 2024.
%
% SEE ALSO:
%   - RSVDinc.m  (incremental randomized SVD)
%   - SRSVDloc.m (blockwise randomized SVD for matrix partitions)
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, October 2022
%   UPC/CIMNE – Barcelona, Spain
%   Email: jhortega@cimne.upc.edu
%--------------------------------------------------------------------------
% Comments by ChatGPT4
% ----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a matrix A, SRSVD returns   a factorization  [U,S,V] of the form
%   A =  U*S*V^T + E   (E = 0 when no truncation is introduced)
% --------------------------
% INPUTS
%------------------------------------------------------------------------------
%  A:    Matrix to be factorized. It can be given as a single numeric array,
%       or as a  1 x q cell array containing a conforming partition of A, i.e.:
%
%     A = {A_1 , A_2 ... A_q}
%
%     In turn, each entry A_i may be a numeric array, or the name of a MAT-file containing
%     the numeric array. The second option is preferred when the whole
%     matrix does not fit into fast memory
%
% -------------------------------------------------------------------------------
%  RELTOL :
%     If  0< RELTOL <1 is a number, RELTOL
%     indicates the error threshold in approximating the whole matrix
%     (in the frobenius norm), i.e.
%     norm(E,'fro')/norm(A,'fro') <= epsilon, where E = A - U*S*V^T
%
%     If RELTOL is a vector then
%     norm(E{i},'fro')/norm(A{i},'fro') <= epsilon(i)
%
% -----------------------------------------------------------------------------
%
%  OUTPUTS
%  ------
%  U --> (M x r)  Matrix of left singular vectors (it approximately spans
%  the column space of A); r denotes the rank of the approximation.
%  V --> (N x r)  Matrix of right singular vectors
%  S ---> Vector of singular values  (r x 1)
%
%  Written by Joaquín A. Hernández Ortega, October  2022
%  UPC/CIMNE. jhortega@cimne.upc.edu

if nargin == 0
    load('tmp3.mat')
end
if nargin <=2
    DATA = [] ;
end
if nargin == 1
    RELTOL = 0 ; 
end

DATA = DefaultField(DATA,'WeightsPremultipy_matrix',[]); % = sqrt(W) ;

if iscell(A)
    if size(A,1) == 1
        if size(A,2) > 1
            [U,S,V] = SRSVDloc(A,RELTOL,DATA);
        else
            %   [U,S,V,eSVD,Rsup] = RSVDT(A,e0,mu,R,DATA)
            DATA.RELATIVE_SVD = 1;
            if ~isempty(DATA.WeightsPremultipy_matrix)
                A{1} = bsxfun(@times,A{1},DATA.WeightsPremultipy_matrix) ;
            end
            
            [U,S,V] = RSVDinc(A{1},RELTOL,[],0,DATA);
        end
    else
        error('Only column block matrices are allowed')
    end
else
    DATA.RELATIVE_SVD = 1;
    if ~isempty(DATA.WeightsPremultipy_matrix)
        A = bsxfun(@times,A,DATA.WeightsPremultipy_matrix) ;
    end
    [U,S,V] = RSVDinc(A,RELTOL,[],0,DATA);
end