function D_QrotINIall = QrotINIall_GLOBAL(QrotINI,INDEXsparseROTmat) 
% Construction of global expandend rotation matrix (diagonal)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx

% Latex notation /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROT.tex
%  \DiagC{\QrotINIall}  =  \diagOL{\QrotINIallE{1}}{\QrotINIallE{2}}{\cdots}{\QrotINIallE{\nelemC}}
%and 
%  \QrotINIallE{e} \defeq \matcdos{\DiagONEc{\QrotINIe{\nelemC}}}{\zero}{\zero}{\ident}
%  INPUT DATA
% MATRIX nSD*nelemC x nSD containing the rotation matrices of the nelemC
% elements 
% INDEXsparseROTmat: Array of indices for expediting the assembly of the
% output, sparse matrix D_QrotINIall. This array is constructed in turn in 
% B_N_matricesEIFEbubCOROT_LRss.m  and QrotINIall_Element.m
% 
% JAHO, 24-oct-2024, thursday, UPC, Campus Nord, Barcelona. 
if nargin == 0
    load('tmp3.mat')
end

% BOUNDARY DOFS
% 
nSD = size(QrotINI,2) ; 
nelem = INDEXsparseROTmat.nelem ; 

% Number of nonzero entries 
%nzmax = length(INDEXsparseROTmat.COLS.DOFsBUB) + length(INDEXsparseROTmat.COLS.DOFsB) ; 
n  = length(INDEXsparseROTmat.COLS.DOFsB)/nSD  + length(INDEXsparseROTmat.COLS.DOFsBUB) ;   % Number of rows/columns
nnodeB = length(INDEXsparseROTmat.COLS.DOFsB)/nSD^2/nelem ; 
%D_QrotINIall = sparse(n,n,nzmax) ; % Allocating space
sVAL_b = zeros(size(INDEXsparseROTmat.COLS.DOFsB))  ; 

for jdim = 1:nSD
    % jdim = Column indexes, local     
    for idim = 1:nSD 
    % idim = Row indexes, local         
    QrotINI_local = repmat(QrotINI(idim:nSD:end,jdim),1,nnodeB) ; % Value of rotation matrices, indices (idim,jdim)    
    QrotINI_local = QrotINI_local' ; 
    indGLO_sparse_ini = (jdim-1)*nSD + idim ; 
    sVAL_b(indGLO_sparse_ini:nSD^2:end) = QrotINI_local(:); 
    end
end


sVAL = [sVAL_b; ones(size(INDEXsparseROTmat.COLS.DOFsBUB(:)))] ; 
iCOL = [INDEXsparseROTmat.COLS.DOFsB; INDEXsparseROTmat.COLS.DOFsBUB(:)] ; 
iROW = [INDEXsparseROTmat.ROWS.DOFsB; INDEXsparseROTmat.ROWS.DOFsBUB(:)] ; 
nzmax = length(sVAL) ; 
D_QrotINIall = sparse(iROW,iCOL,sVAL,n,n,nzmax) ; 
