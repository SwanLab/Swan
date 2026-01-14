function [ iROWS,jCOLS] = FintCunassemb_Element(e,number_accumulated_DOFs,EIFEoper_all)
%--------------------------------------------------------------------------
%  FintCunassemb_Element
%
%  Computes the row and column index vectors required to assemble the global
%  extended internal force vector in a sparse format for a given element `e`.
%
%  This corresponds to the construction of the global operator:
%
%      \DiagC{\FintCunassemb} = [Fint^1; Fint^2; ...; Fint^nElem]
%
%  where each elemental contribution (composed of boundary and bubble DOFs)
%  is stacked vertically into a global column-wise structure, and
%  assembly is carried out using `iROWS` and `jCOLS`.
%
%  INPUTS:
%    - e                    : Element number (1-based)
%    - number_accumulated_DOFs : Offset from previous elements (used to place
%                                this element’s DOFs in global matrix)
%    - EIFEoper_all         : Structure containing:
%         > INFO.DOFsB    : Local indices of boundary DOFs
%         > INFO.DOFsBUB  : Local indices of bubble DOFs
%
%  OUTPUTS:
%    - iROWS : Global row indices for sparse matrix assembly
%    - jCOLS : Global column indices (set to element number `e`)
%
%  The pair (iROWS, jCOLS) is later used in the global block structure
%  of the extended internal force vector and associated operators
%  such as:
%      - D_FintCunassemb_Element
%      - D_BmatIst * d  or  QrotINIallEᵗ * f_int^e
%
%  USAGE:
%    Typically used in the function `B_N_matricesEIFEbubCOROT_LRss.m`
%    and other EIFEM precomputations related to global operator assembly.
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 6-Feb-2025 (Thursday)
%    Comments by ChatGPT4, 13-May-2025
%
%--------------------------------------------------------------------------



% Determination of indices for constructing for D_FintCunassemb_Element 
% Latex code:   \DiagC{\FintCunassemb} 
% iROWS  --- Indices Rows
% jCOLS --- Indices Columns 
% JAHO, 6-feb-2025, UPC, BARCELONA, Thursday
if nargin == 0
    load('tmp1.mat')
end

% BOUNDARY DOFs/BUBBLE DOFs
DOFsB = EIFEoper_all.INFO.DOFsB ; % Boundary DOFs, local index
DOFsBUB = EIFEoper_all.INFO.DOFsBUB ; % Bubble DOFs, local indexing

DOFsLOC = [DOFsB,DOFsBUB]; 

% For e = 1 
% irows = DOFsLOC  ;   jCOLS = e*ones(size(DOFsLOC)) ; 
% For e 
% irows = DOFsLOC + number_accumulated_DOFs ;    jCOLS = e*ones(size(DOFsLOC)) ;
 
iROWS = DOFsLOC(:) + number_accumulated_DOFs  ; 
jCOLS =  e*ones(size(DOFsLOC))   ;
jCOLS = jCOLS(:); 


% nSD = 1; % Number of spatial dimension
% % For 2D, the excerpt below simply sets:  iROWS = [1,2,1,2], jCOLS = [1,1,2,2] ;
% iROWS = repmat([1:nSD],1,nSD) ;
% jCOLS = reshape(iROWS,nSD,nSD)' ;
% jCOLS = jCOLS(:)' ;
% %- Now we make nnodesB copies of iROWS and jCOLS
% iROWSallE = repmat(iROWS,nnodesB,1) ;
% jCOLSallE = repmat(jCOLS,nnodesB,1) ;
% % Each of these rows have to be summed up sequentially by 0, 2,4... ...
% % nnodesB
% iROWSallE = bsxfun(@plus,iROWSallE,[0:nSD:(nnodesB-1)*nSD]')' ;
% iROWSallE_b = iROWSallE(:) ;
% jCOLSallE = bsxfun(@plus,jCOLSallE,[0:nSD:(nnodesB-1)*nSD]')' ;
% jCOLSallE_b = jCOLSallE(:) ;
% % Finally, the matrices in a sparse format (just one colum)
% sMAT_b = repmat(QrotINI(:),nnodesB,1) ;
% % Now we need to add the bubble DOFs. This is simply a identity matrix
% %
% %  Thus
% iROWSallE = [iROWSallE_b;DOFsBUB] ; 
% jCOLSallE = [jCOLSallE_b;DOFsBUB] ; 
% sMAT = [sMAT_b; ones(length(DOFsBUB))'] ; 
% iROWS = [] ; 
% iROWS.DOFsB = iROWSallE_b ; 
% iROWS.DOFsBUB = DOFsBUB ; 
% jCOLS = [] ; 
% jCOLS.DOFsB = jCOLSallE_b ; 
% jCOLS.DOFsBUB = DOFsBUB ; 
% 
% QrotINIallE = sparse(iROWSallE,jCOLSallE,sMAT) ;
