function   [i,j] = IndicesCtang(m,p)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: IndicesCtang
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   Generates the row and column indices `(i, j)` for assembling the consistent tangent stiffness
%   matrix (`Ctang`) in a block-wise manner. The resulting indices are used to populate sparse 
%   matrices corresponding to Gauss-point-level tangent operators.
%
%   In particular, the function returns the global indices of local (element-level) contributions 
%   for multiple Gauss points, assuming that each Gauss point contributes a `p × p` block.
%
% USAGE:
%   [i, j] = IndicesCtang(m, p)
%
% INPUTS:
%   - m : Total number of rows/columns in the final sparse matrix (e.g., `m = ngaus * p`)
%   - p : Size of each block (typically the number of strain components)
%
% OUTPUTS:
%   - i : Row indices of the non-zero entries in the sparse matrix
%   - j : Column indices corresponding to `i`
%
% FUNCTIONALITY:
%   - Constructs all index pairs `(ielem, jelem)` for a `p × p` matrix.
%   - Repeats this index pattern `ngaus = m/p` times, offsetting each by a block.
%   - Useful when assembling tangent operators from Gauss-point-based contributions:
%       Ctang(i(k), j(k)) = Ctang(i(k), j(k)) + local_entry(k)
%
% EXAMPLE:
%   For `p = 3`, `m = 6`:
%     - `ngaus = 2`, number of Gauss points
%     - The function returns indices for two 3×3 blocks at positions 1:3 and 4:6
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE  
%   Date: Unknown (template code structure), revised preamble May-2025
%   Comments by ChatGPT4, 12-May-2024
% DEPENDENCIES:
%   - None (uses built-in MATLAB functions only)
%
% ---------------------------------------------------------------------------------------------------


if nargin == 0
    p = 2;
    m = 6 ;
end

n = m ;
ngaus = m/p ;
nzmax = m*p ;
ij = zeros(p^2,2) ;
for ielem = 1:p
    for jelem = 1:p
        iniELEM = (ielem-1)*p +jelem;
        ij(iniELEM,1) = ielem ;
        ij(iniELEM,2) = jelem ;
    end
end

ij = repmat(ij,ngaus,1);
indSUM = 0:p:(ngaus-1)*p ;
indSUM = repmat(indSUM,p^2,1);
indSUM = reshape(indSUM,size(indSUM,1)*size(indSUM,2),1) ;
ij = bsxfun(@plus,ij,indSUM) ;
i = ij(:,1);
j = ij(:,2);