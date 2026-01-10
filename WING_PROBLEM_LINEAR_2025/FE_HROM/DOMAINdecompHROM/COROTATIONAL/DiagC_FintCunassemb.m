function D_FintCunassemball = DiagC_FintCunassemb(FintCunassemb,INDEXsparseFINT)  
% Construction of  
%  \DiagC{\FintCunassemb} \defeq \diagOL{\FintCe{1}}{\FintCe{2}}{\cdots}{\FintCe{\nelemC}}
% % See  Small strains/Large rotations (or Small rotations/Large strains)
% JAHO, 15-feb-2025, SATURDAY, 8:18, Balmes 185,  Barcelona.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% and
%/home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf

if nargin == 0
    load('tmp3.mat')
end

% BOUNDARY DOFS
% 
nSD = 1 ; 
nelem = unique(INDEXsparseFINT.COLS) ; 
nrows= length(INDEXsparseFINT.COLS) ; 
ncols = length(nelem) ; 
nzmax = nrows ; 
D_FintCunassemball = sparse(INDEXsparseFINT.ROWS,INDEXsparseFINT.COLS,FintCunassemb,nrows,ncols,nzmax) ; 


% 
% % Number of nonzero entries 
% %nzmax = length(INDEXsparseROTmat.COLS.DOFsBUB) + length(INDEXsparseROTmat.COLS.DOFsB) ; 
% n  = length(INDEXsparseROTmat.COLS.DOFsB)/nSD  + length(INDEXsparseROTmat.COLS.DOFsBUB) ;   % Number of rows/columns
% nnodeB = length(INDEXsparseROTmat.COLS.DOFsB)/nSD^2/nelem ; 
% %D_FintCunassemball = sparse(n,n,nzmax) ; % Allocating space
% sVAL_b = zeros(size(INDEXsparseROTmat.COLS.DOFsB))  ; 
% 
% for jdim = 1:nSD
%     % jdim = Column indexes, local     
%     for idim = 1:nSD 
%     % idim = Row indexes, local         
%     FintCunassemb_local = repmat(FintCunassemb(idim:nSD:end,jdim),1,nnodeB) ; % Value of rotation matrices, indices (idim,jdim)    
%     FintCunassemb_local = FintCunassemb_local' ; 
%     indGLO_sparse_ini = (jdim-1)*nSD + idim ; 
%     sVAL_b(indGLO_sparse_ini:nSD^2:end) = FintCunassemb_local(:); 
%     end
% end
% 
% 
% sVAL = [sVAL_b; ones(size(INDEXsparseROTmat.COLS.DOFsBUB))] ; 
% iCOL = [INDEXsparseROTmat.COLS.DOFsB; INDEXsparseROTmat.COLS.DOFsBUB] ; 
% iROW = [INDEXsparseROTmat.ROWS.DOFsB; INDEXsparseROTmat.ROWS.DOFsBUB] ; 
% nzmax = length(sVAL) ; 
