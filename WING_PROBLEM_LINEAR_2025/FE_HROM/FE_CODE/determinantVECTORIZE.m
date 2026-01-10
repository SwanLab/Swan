function    detJ= determinantVECTORIZE(Jinp,ndim)
%%%%
% -------------------------------------------------------------------------
% COMMENTS (generated automatically by ChatGPT on 7-Nov-2025)
%
% PURPOSE:
%   Vectorized determinant of stacked Jacobian blocks.
%   Given a tall matrix
%       Jinp = [J_1; J_2; ...; J_nelem],
%   where each J_e is an (ndim x ndim) block (arranged one after another in
%   rows), this routine returns detJ, a column vector with det(J_e) for all
%   elements e = 1..nelem.
%
% INPUTS:
%   Jinp : ((nelem*ndim) x ndim) matrix formed by vertically stacking the
%          per-element Jacobians J_e. The e-th block occupies rows
%          ((e-1)*ndim+1) : (e*ndim).
%   ndim : spatial dimension and block size. Supported values: 1, 2, or 3.
%
% OUTPUTS:
%   detJ : (nelem x 1) vector of determinants det(J_e).
%
% DETAILS / IMPLEMENTATION NOTES:
%   - For ndim == 1, the determinant is the scalar itself → detJ = Jinp.
%   - For ndim == 2 or 3, the routine:
%       1) Re-slices Jinp into a small cell array J{i,j}, where each cell is
%          a (nelem x 1) column collecting the (i,j) entry across all blocks
%          (i.e., rows i, i+ndim, i+2*ndim, ... of column j in Jinp).
%       2) Applies the explicit closed-form determinant formula componentwise:
%            * 2D: det = a11*a22 - a12*a21
%            * 3D: det = a11*a22*a33 - a11*a23*a32 - a12*a21*a33
%                        + a12*a23*a31 + a13*a21*a32 - a13*a22*a31
%   - This avoids per-element loops and uses elementwise array operations.
%
% ASSUMPTIONS & SANITY:
%   - Jinp must be correctly stacked (no padding rows).
%   - No reorientation checks are done here; negative detJ values are
%     returned as-is (caller may warn/error if orientation is invalid).
%
% EXAMPLE (self-test path when nargin==0):
%   Builds three 3x3 blocks J1, 4*J1, 6*J1, stacks them in Jinp, and
%   compares the output detJ with det(Jk).
%
% AUTHOR / HISTORY:
%   Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Oct-2015
%   Comments clarification (this header): 7-Nov-2025
% -------------------------------------------------------------------------

% Given a matrix Jinp = [Jinp_1; Jinp_22 ... ; Jinp_nelem]  (Jinp_i is ndim x ndim),
% determinantVECTORIZE returns a vector consisting of the determinants of
% each block matrix Jinp_i
% ndim = 2,3
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Oct-2015
if nargin == 0
    J1 = [3 5 6; 7 8 7 ; 0 0 3] ;  d1 =det(J1) ;
    J2 = 4*[3 5 6; 7 8 7 ; 0 0 3] ;  d2 =det(J2) ;
    J3 = 6*[3 5 6; 7 8 7 ; 0 0 3] ;  d3 =det(J3) ;
    Jinp = [J1; J2; J3];
    ndim = 3 ;
end
%

if   ndim == 1
    detJ = Jinp ;
else
    
    J = cell(ndim,ndim) ;
    for i=1:ndim
        iglo = i:ndim:size(Jinp,1) ;
        for j=1:ndim
            jglo = j ;
            J{i,j} = Jinp(iglo,jglo) ;
        end
    end
    if  ndim==2
        detJ = J{1,1}.*J{2,2}  -   J{1,2}.*J{2,1} ;
    elseif ndim == 3
        detJ = J{1,1}.*J{2,2}.*J{3,3} - J{1,1}.*J{2,3}.*J{3,2} - J{1,2}.*J{2,1}.*J{3,3} + J{1,2}.*J{2,3}.*J{3,1} ...
            + J{1,3}.*J{2,1}.*J{3,2} - J{1,3}.*J{2,2}.*J{3,1} ;
        
    end
    
end
