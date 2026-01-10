function Ahinv = MultiUpdateInverseHermitian(Binv,jrowMAT)
%--------------------------------------------------------------------------
% function Ahinv = MultiUpdateInverseHermitian(Binv, jrowMAT)
%
% PURPOSE:
%   Performs a **recursive update of the inverse** of a Hermitian matrix `Binv`
%   by sequentially eliminating the rows and columns specified in `jrowMAT`.
%   It is based on the Sherman–Morrison–Woodbury update for Hermitian matrices.
%
%   This is particularly useful in the context of **basis compression**, 
%   hyperreduction techniques (e.g., DEIM, ECM), or adaptive algorithms 
%   where certain basis functions or columns are removed to improve efficiency.
%
% INPUTS:
%   - Binv     : Initial inverse matrix of a Hermitian matrix B = AᵗA ∈ ℝⁿˣⁿ
%   - jrowMAT  : Vector of indices (row/col positions) to be recursively removed
%
% OUTPUTS:
%   - Ahinv    : Updated inverse matrix of reduced system (with rows/cols removed)
%
% INTERNAL PROCEDURE:
%   - Sorts `jrowMAT` in ascending order to avoid index mismatches.
%   - Applies `UpdateInverseHermitian` for each index `jrow`, corrected as:
%       jrow = jrowMAT(i) - i + 1
%     since after each removal, the matrix shrinks.
%   - The update is applied recursively.
%
% EXAMPLE USAGE (from test case in nargin == 0 block):
%   m = 10; n = 7;
%   jrowMAT = [8 9];
%   A = randn(m,n); a = randn(m,length(jrowMAT));
%   Bor = zeros(m,n+length(jrowMAT));
%   Bor(:,setdiff(1:end,jrowMAT)) = A; Bor(:,jrowMAT) = a;
%   Binv = inv(Bor'*Bor);
%   Ahinv = MultiUpdateInverseHermitian(Binv, jrowMAT);
%
% DEPENDENCIES:
%   - Requires function `UpdateInverseHermitian` (not shown here) to perform
%     the single-column/row update of Hermitian matrix inverses.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, 24 March 2017.
%--------------------------------------------------------------------------

% Recursive application of UpdateInverseHermitian
%  J.A. Hdez, 24 March 2017
if nargin == 0
    m = 10; n =7; 
    jrowMAT =[8 9] ; 
    A = randn(m,n) ; a = randn(m,length(jrowMAT)) ;
    Bor = zeros(m,n+length(jrowMAT)) ; 
    indOR = 1:size(Bor,2) ; 
    indOR(jrowMAT) = [] ; 
    Bor(:,indOR) = A ; 
    Bor(:,jrowMAT) = a;     
    B = [Bor'*Bor] ;  
    Binv = inv(B);     
    AhinvREAL = inv(A'*A)      ;
        addpath('/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/FE_CODE/')

end
 
jrowMAT = sort(jrowMAT) ; 
BinvOLD = Binv ; 
for i = 1:length(jrowMAT)
    jrow = jrowMAT(i)-i+1 ; 
    Ahinv = UpdateInverseHermitian(BinvOLD,jrow) ; 
    BinvOLD = Ahinv ; 
end

