function [NDOFS_pernode,MESHextended_COOR,MESHextended_CN,...
    DOFS_TO_KEEP,DOFS_bLOC,DOFS_bubLOC,DOFS_ghost] = GhostDOFs(COOR,CN,DATA,nDOFSbub)
%--------------------------------------------------------------------------
%  GhostDOFs
%
%  Constructs the **extended mesh structure** required to account for bubble modes
%  in the EIFEM formulation. It returns updated coordinates and connectivities that include
%  pseudo-centroid (bubble) nodes, and classifies all DOFs into boundary, bubble, and ghost.
%
%  INPUTS:
%    - COOR         : coordinates of original (boundary) nodes
%    - CN           : connectivity matrix of coarse-scale elements (only boundary nodes)
%    - DATA         : structure for internal settings (may be enriched inside)
%    - nDOFSbub     : number of DOFs associated to bubble nodes (per node)
%
%  OUTPUTS:
%    - NDOFS_pernode    : number of DOFs per node (including extended ones)
%    - MESHextended_COOR: extended coordinate matrix including bubble nodes
%    - MESHextended_CN  : extended connectivity matrix
%    - DOFS_TO_KEEP     : vector of indices for DOFs to keep (i.e., excluding ghost DOFs)
%    - DOFS_bLOC        : local indices of retained boundary DOFs
%    - DOFS_bubLOC      : local indices of retained bubble DOFs
%    - DOFS_ghost       : DOFs to be eliminated during assembly (e.g., due to zero displacement)
%
%  THEORY:
%    - Each coarse-scale element is enriched with one pseudo-centroid node (bubble node)
%      whose coordinates are the average of its boundary node coordinates.
%    - The corresponding DOFs are allocated at the end of the coordinate array, and are
%      distinguished from boundary DOFs for assembly purposes.
%    - Depending on the number of DOFs assigned to bubble nodes versus spatial dimension:
%        • If nDOFSbub > ndim: excess DOFs in boundary nodes are marked as ghost
%        • If ndim > nDOFSbub: excess DOFs in bubble nodes are marked as ghost
%        • If equal: no ghost DOFs are present
%
%  ROLE IN EIFEM:
%    - This construction enables the separation of bubble and boundary DOFs for:
%        • Tangent matrix assembly (e.g., block partitions)
%        • Filtering of non-physical variables in reduced-order expansions
%    - Ghost DOFs are purged from the final system using `DOFS_TO_KEEP`
%    - The extended coordinate matrix is used in building the matrix `XcALL` and in
%      transformations such as:
%         Qᵀ (dc - d̂cbQ)   (Eq. 554 in *EIFEM_largeROTfinal.pdf*):contentReference[oaicite:0]{index=0}
%
%  REFERENCES:
%    - Section 12.6, Eq. (554), Eq. (655), and Appendix 17 of *EIFEM_largeROTfinal.pdf*
%    - Implementation concept from `BubbleModeImplement.mlx`
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 26–28 Oct 2024
%    Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------





% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/BubbleModeImplement.mlx

%  EXTENDED COORDINATE MATRIX
%  --------------------------
% PSEUDO-CENTROIDS OF THE ELEMENT
% -------------------------------
ndimSP = size(COOR,2) ;  % Number of spatial dimensions

if nDOFSbub == 0
    NDOFS_pernode = ndimSP ; 
    MESHextended_COOR = COOR; 
    MESHextended_CN = CN ; 
    DOFS_TO_KEEP = 1:prod(size(COOR)); 
    DOFS_bLOC = DOFS_TO_KEEP ; 
    DOFS_bubLOC = [] ; 
    DOFS_ghost = [] ; 
    
else

COORcent = zeros(size(CN,1),size(COOR,2)) ;
for inodeE = 1:size(CN,2)
    NODESloc = CN(:,inodeE) ;
    COORcent = COORcent + COOR(NODESloc,:) ;
end
COORcent  = COORcent/size(CN,2)   ;
% Extended coordinate matrix / connectivities

columCENTROID_CN = size(COOR,1) + (1:size(CN,1))' ;

MESHextended_COOR = [COOR; COORcent] ;
MESHextended_CN = [CN,columCENTROID_CN] ;

% DOFs to be eliminated (GHOST_DOFS)
% -----------------------------------------------------------------

if nDOFSbub > ndimSP
    DATA.MESHextended.NDOFS_pernode = nDOFSbub ;
    % Extended set of DOFs
    NDOFS_pernode =  nDOFSbub;
    NDOFS_TOTAL = size(MESHextended_COOR,1)*NDOFS_pernode ;
    NDOFS_bnd = size(COOR,1)*NDOFS_pernode ;
    
    DOFS_bnd = 1:NDOFS_bnd ;
    DOFS_bnd = reshape(DOFS_bnd,NDOFS_pernode,[]) ;
    
    DOFS_ghost = DOFS_bnd(ndimSP+1:end,:) ;
    DOFS_b  = DOFS_bnd(1:ndimSP ,:) ;
    DOFS_bub = ((size(COOR,1)*NDOFS_pernode)+1):NDOFS_TOTAL ;
    
    %DATA.MESHextended.DOFS_bub = DOFS_bub(:) ;
    %DATA.MESHextended.DOFS_b = DOFS_b(:) ;
    DOFS_ghost = DOFS_ghost(:) ;
    
    % thus
    
    
elseif ndimSP > nDOFSbub
    DATA.MESHextended.NDOFS_pernode = ndimSP ;
    % Extended set of DOFs
    NDOFS_pernode =  ndimSP;
    NDOFS_TOTAL = size(MESHextended_COOR,1)*NDOFS_pernode ;
    NDOFS_bnd = size(COOR,1)*NDOFS_pernode ;
    
    DOFS_bnd = 1:NDOFS_bnd ;
    DOFS_b = DOFS_bnd(:) ;
    
    DOFS_bub = ((size(COOR,1)*NDOFS_pernode)+1):NDOFS_TOTAL ;
    DOFS_bub = reshape(DOFS_bub,NDOFS_pernode,[]) ;
    DOFS_ghost = DOFS_bub(nDOFSbub+1:end,:) ;
    DOFS_bub  = DOFS_bub(1:nDOFSbub ,:) ;
    
    %     DATA.MESHextended.DOFS_bub = DOFS_bub(:) ;
    %     DATA.MESHextended.DOFS_b = DOFS_b(:) ;
     DOFS_ghost = DOFS_ghost(:) ;
    
    
else
    DATA.MESHextended.NDOFS_pernode = ndimSP ;
    NDOFS_pernode =  ndimSP;
    NDOFS_TOTAL = size(MESHextended_COOR,1)*NDOFS_pernode ;
    NDOFS_bnd = size(COOR,1)*NDOFS_pernode ;
    
    DOFS_bnd = 1:NDOFS_bnd ;
    DOFS_b = DOFS_bnd(:) ;
    
    DOFS_bub = ((size(COOR,1)*NDOFS_pernode)+1):NDOFS_TOTAL ;
    
    %     DATA.MESHextended.DOFS_bub = DOFS_bub(:) ;
    %     DATA.MESHextended.DOFS_b = DOFS_b(:) ;
     DOFS_ghost = [];
    
end

DOFS_all= 1:NDOFS_TOTAL ;
[DOFS_TO_KEEP,AAA ]= setdiff(DOFS_all(:), DOFS_ghost) ;

[dummy1, DOFS_bLOC]  = intersect( DOFS_TO_KEEP(:),DOFS_b(:)) ;

[dummy2, DOFS_bubLOC]  = intersect( DOFS_TO_KEEP(:),DOFS_bub(:)) ;

end

% % So, for instance, in the case used for testing this routine 
% %
% DOFS_bLOC'
% 
% ans =
% 
%      1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20...
%       119   120   121   122   123   124   125   126   127   128
% and 
%   DOFS_bubLOC'
% 
% ans =
% 
%    129   130   131   132   133   134   135   136   137   138    .....
% 
%    318   319   320   321   322   323   324

% But this might not be always the case (I mean, they might not be consecutive)...For instance, it won't be like this
% for quadratic elements 