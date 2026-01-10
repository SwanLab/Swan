function [LboolCall,XcALL] =   BooleanOperatorCoarse(COOR,CN,number_accumulated_DOFs,NumberBubbleDOFS_perELEMENT)
%--------------------------------------------------------------------------
%  BooleanOperatorCoarse
%
%  Constructs the global **Boolean matrix** `LboolCall` for the coarse-scale
%  displacement vector in EIFEM, and computes the block-wise coordinate vector
%  `XcALL` used in the rotation update and geometric terms.
%
%  This function returns the global matrix:
%
%     LboolCall := [L_c(1); L_c(2); ...; L_c(nelem)]
%
%  where each `L_c(e)` extracts the local vector of coarse-scale DOFs (interface + bubble)
%  from the global displacement vector `dc`, as described in:
%
%     d_c(e) = L_c(e) * d_c    (Eq. 446 and 535 in *EIFEM_largeROTfinal.pdf*):contentReference[oaicite:0]{index=0}
%
%  INPUTS:
%    - COOR                          : Coordinates of boundary and bubble nodes
%    - CN                            : Connectivity matrix of boundary nodes
%    - number_accumulated_DOFs       : Total number of DOFs across all elements
%    - NumberBubbleDOFS_perELEMENT   : Vector with the number of bubble DOFs per element
%
%  OUTPUTS:
%    - LboolCall : Global Boolean matrix for coarse-scale displacement assembly
%    - XcALL     : Vector of coarse-scale nodal coordinates (augmented with zeros for bubble DOFs)
%
%  NOTES:
%    - The Boolean matrix structure reflects the block decomposition:
%
%         [d_cb^e ; d_cbub^e] = [L_cb^e  0; 0  L_cbub^e] * [d_cb ; d_cbub]
%
%      which is consistent with Eq. (538) in *EIFEM_largeROTfinal.pdf*:contentReference[oaicite:1]{index=1}
%
%    - The bubble DOFs are indexed globally after the boundary DOFs, and are mapped
%      using local accumulation offsets (`ndofACUM_bub`).
%
%    - The vector `XcALL` is also built element-wise, stacking boundary node coordinates
%      (flattened column-wise) and appending zeros for bubble DOFs. It is used in
%      expressions such as:
%
%         d̂c′_Q = XcALLloc - Qᵀ XcALL   (Eq. 554):contentReference[oaicite:2]{index=2}
%
%  ROLE:
%    - `LboolCall` is used during the assembly of the internal force vector:
%         Fint_c = LboolCallᵀ * Fint_local   (Eq. 534):contentReference[oaicite:3]{index=3}
%
%    - This operator does not change during Newton iterations and can be precomputed.
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Versions: 26-Oct-2024 / 28-Oct-2024
%    Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------



% --------------------------------------------------------------------
% EXTENDED CONNECTIVITY MATRIX, BOOLEAN OPERATOR 
% --------------------------------------------------
% We know at this point that the   number of DOFs is  number_accumulated_DOFs
% We also know the number of bubble DOFs per element  NumberBubbleDOFS_perELEMENT
% Our strategy will go as follows. Suppose we define a coordinate matrix in
% which the bubble nodes are placed at the end; 
% COOR = [COOR_boundary; COOR_bubble], so that 
% COOR_bubble(ielem) is the bubble node of the ielem-th element 
%  This function returns the global Boolean matrix defined by (see
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx)
%  /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROT.tex

%  \LboolCall \defeq  \colcuatro{\LboolC{1}}{\LboolC{2}}{\vdots}{\LboolC{\nelemC}}
% where
%  \coldos{\dCbE{e}}{\dCbubE{e}} = \matcdos{\LboolCb{e}}{\zero}{\zero}{\LboolCbub{e}} \coldos{\dCb}{\dCbub}

% Coordinate vector is also computed here (XcALL)

% JAHO, 26-Oct-2024, Saturday, Balmes 185, Barcelona 
%       28-Oct-2024,   Monday, Balmes 185, Barcelona
% ---------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

nnodeBe = size(CN,2) ;  % Number of interface boundary nodes per element
nSD = size(COOR,2) ; 
nnodeB = size(COOR,1); % Total number of boundary nodes  
nDOFSc = number_accumulated_DOFs ;  % Total number of (coarse-scale ) DOFs  (nnodeB*nSD + size(CN,1))
nDOFSbub = sum(NumberBubbleDOFS_perELEMENT) ; 
nelem = size(CN,1) ; 
nDOFSb_e = nnodeBe*nSD ; % Number of boundary DOFs per element 
nDOFSb = nnodeB*nSD ;    % Total number of boundary DOFs 
LboolCall = cell(nelem,1) ; % Boolean matrix (coarse-scale assembly operator)
XcALL = cell(nelem,1) ; % Coordinate vector   

ndofACUM_bub = 0 ; 
for ielem = 1:nelem 
      % --------------------------------
      % Boolean matrix boundary DOFs
      % --------------------------------
      nodesGLO  = CN(ielem,:) ;  
      COORnodes = COOR(nodesGLO,:) ; 
      COORnodes = COORnodes' ; 
      COORnodes = COORnodes(:) ; 
      
      
      
      ndofsGLO = small2large(nodesGLO,nSD) ; 
      LboolCb = sparse(1:nDOFSb_e,ndofsGLO,ones(size(ndofsGLO)),nnodeBe*nSD,nnodeB*nSD) ;
      % --------------------------------
      % Boolean matrix bubble DOFs
      % --------------------------------
      ndofsBUBelem = NumberBubbleDOFS_perELEMENT(ielem) ; 
      idofsGLO = (ndofACUM_bub+1):(ndofACUM_bub+ndofsBUBelem) ;  
      LboolCbub = sparse(1:ndofsBUBelem,idofsGLO,ones(size(idofsGLO)),ndofsBUBelem,nDOFSbub) ;
      
      LboolCall{ielem} = blkdiag(LboolCb,LboolCbub) ; 
      ndofACUM_bub = ndofACUM_bub + ndofsBUBelem ; 
      
      XcALL{ielem} = [COORnodes; zeros(ndofsBUBelem,1)] ; 
end 


LboolCall = cell2mat(LboolCall) ; 
XcALL = cell2mat(XcALL); 