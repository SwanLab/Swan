function [U,S,V] = RSVDTrowblock(A,RELTOL,DATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a matrix A, RSVDTrowblock returns   a factorization  [U,S,V] of the form
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
    load('tmp.mat')
end
if nargin == 2
    DATA = [] ;
end

DATA = DefaultField(DATA,'WeightsPremultipy_matrix',[]); % = sqrt(W) ;

if iscell(A)
    if size(A,1) == 1
        if size(A,2) > 1
            [U,S,V] = RSVDTrowblockLOC(A,RELTOL,DATA);
        else
            %   [U,S,V,eSVD,Rsup] = RSVDT(A,e0,mu,R,DATA)
            DATA.RELATIVE_SVD = 1;
            if ~isempty(DATA.WeightsPremultipy_matrix)
                A{1} = bsxfun(@times,A{1},DATA.WeightsPremultipy_matrix) ;
            end
            
            [U,S,V] = RSVDT(A{1},RELTOL,[],0,DATA);
        end
    else
        error('Only row block matrices are allowed')
    end
else
    DATA.RELATIVE_SVD = 1;
    if ~isempty(DATA.WeightsPremultipy_matrix)
        A = bsxfun(@times,A,DATA.WeightsPremultipy_matrix) ;
    end
    [U,S,V] = RSVDT(A,RELTOL,[],0,DATA);
end