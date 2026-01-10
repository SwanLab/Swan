function [QrotINIallE,iROWS,jCOLS] = QrotINIall_Element(EIFEoper_all,QrotINI)
%--------------------------------------------------------------------------
%  QrotINIall_Element
%
%  Constructs the block-diagonal **initial rotation matrix** for an EIFEM
%  element by assembling the tensor product of the element’s initial rotation
%  matrix (`QrotINI`) across all **boundary DOFs**. Bubble DOFs are assumed
%  to remain invariant under rotation and are thus assigned an identity matrix.
%
%  This function forms the matrix:
%
%      QrotINIallE{e} := [ kron(I_nodesB, QrotINI)    0 ;
%                           0                         I_bubble ]
%
%  which is used to rotate element-level vectors and matrices into a common
%  reference frame (see Eq. (32) in *EIFEM_largeROTfinal.pdf*).
%
%  INPUTS:
%    - EIFEoper_all : structure containing fields:
%         > INFO.DOFsB    : local indices of boundary DOFs
%         > INFO.DOFsBUB  : local indices of bubble DOFs
%    - QrotINI       : initial rotation matrix of the element (ndim × ndim)
%
%  OUTPUTS:
%    - QrotINIallE : sparse matrix (nDOFs × nDOFs) to be used in element-wise
%                    transformations (e.g., internal force computation)
%    - iROWS, jCOLS : sparse matrix row/column indices for efficient assembly
%       > iROWS.DOFsB    and jCOLS.DOFsB    for boundary DOFs
%       > iROWS.DOFsBUB  and jCOLS.DOFsBUB  for bubble DOFs (identity)
%
%  Notes:
%    - This matrix is assembled elementwise and later integrated into
%      the global operator using `INDEXsparseROTmat`.
%    - For each boundary node, QrotINI is replicated blockwise using Kronecker-like structure.
%    - Bubble DOFs are handled independently and contribute only identity.
%
%  Use in:
%    - `B_N_matricesEIFEbubCOROT_LRss.m`
%    - Co-rotational downscaling and upscaling operations
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 24-Oct-2024
%    Comments by CHATGPT4, 13-May-2025
%--------------------------------------------------------------------------




% Determination of initial rotation matrix,   repeated  for all boundary
% DOFs, in EIF elements. Bubble DOFs are not transformed (IDENTITY MATRIX)
% Latex code:  \QrotINIallE{e} \defeq \matcdos{\DiagONEc{\QrotINIe{\nelemC}}}{\zero}{\zero}{\ident}
% iROWS.DOFsB  --- Indices Rows, for sparse assembly
% jCOLS.DOFsB --- Indices Columns (only boundary DOFs)
% JAHO, 24-OCT-2024, UPC, BARCELONA

% BOUNDARY DOFs/BUBBLE DOFs
DOFsB = EIFEoper_all.INFO.DOFsB ; % Boundary DOFs, local index
DOFsBUB = EIFEoper_all.INFO.DOFsBUB ; % Bubble DOFs, local indexing
nSD = size(QrotINI,1) ; % Number of spatial dimension
nnodesB = length(DOFsB)/nSD ; % Number of nodes
% For 2D, the excerpt below simply sets:  iROWS = [1,2,1,2], jCOLS = [1,1,2,2] ;
iROWS = repmat([1:nSD],1,nSD) ;
jCOLS = reshape(iROWS,nSD,nSD)' ;
jCOLS = jCOLS(:)' ;
%- Now we make nnodesB copies of iROWS and jCOLS
iROWSallE = repmat(iROWS,nnodesB,1) ;
jCOLSallE = repmat(jCOLS,nnodesB,1) ;
% Each of these rows have to be summed up sequentially by 0, 2,4... ...
% nnodesB
iROWSallE = bsxfun(@plus,iROWSallE,[0:nSD:(nnodesB-1)*nSD]')' ;
iROWSallE_b = iROWSallE(:) ;
jCOLSallE = bsxfun(@plus,jCOLSallE,[0:nSD:(nnodesB-1)*nSD]')' ;
jCOLSallE_b = jCOLSallE(:) ;
% Finally, the matrices in a sparse format (just one colum)
sMAT_b = repmat(QrotINI(:),nnodesB,1) ;
% Now we need to add the bubble DOFs. This is simply a identity matrix
%
%  Thus
iROWSallE = [iROWSallE_b;DOFsBUB(:)] ; 
jCOLSallE = [jCOLSallE_b;DOFsBUB(:)] ; 
sMAT = [sMAT_b; ones(length(DOFsBUB),1)] ; 
iROWS = [] ; 
iROWS.DOFsB = iROWSallE_b ; 
iROWS.DOFsBUB = DOFsBUB ; 
jCOLS = [] ; 
jCOLS.DOFsB = jCOLSallE_b ; 
jCOLS.DOFsBUB = DOFsBUB ; 

QrotINIallE = sparse(iROWSallE,jCOLSallE,sMAT) ;
