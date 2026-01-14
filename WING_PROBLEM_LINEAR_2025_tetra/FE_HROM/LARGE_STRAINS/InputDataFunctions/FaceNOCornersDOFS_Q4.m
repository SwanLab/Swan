function [b,f,F,alphaB,betaB,Z]  =FaceNOCornersDOFS_Q4(rnodLOC,ndim,MESH,DATALOC,COMB_FACES,GEOproperties)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/13_Q4_nocorners.mlx
% 


if nargin == 0
    load('tmp.mat')
end

% STEP 1) IDENTIFY THE DISTINCT SETS OF THE DOFS  (CASE WITH CORNERS ! )
%  \item To impose periodic boundary conditions, it is necessary to divide the nodes of the interface boundary of the rectangular element in the following subsets
%
%  \begin{equation}
%   \b =  \f_1 \union  \f_2  \union \f_3  \union \f_4    
%   \end{equation}
% 
%   Likewise, we denote the counterpart of the above sets in local numbering by 
%   
%    \begin{equation}
%   \{1,2 \ldots N\} =  \F_1 \union  \F_2  \union \F_3  \union \F_4     
%   \end{equation}
b = small2large(rnodLOC,ndim) ;  % THESE ARE THE BOUNDARY DOFS
FACES_match = DATALOC.PERIODICITY_FACES  ;
   

% FACES F1,F2,F3  and F4
% 
% % subSTEP 1. Remove corners
% % -------------------------
 nfaces =4;
 
% STEPs . Sort nodes FACES such that NODE_1(i) is the opposite of NODE_3(i)



f = cell(nfaces,1) ; 
F = cell(nfaces,1) ; 

for imatch = 1:length(FACES_match)
    iface = FACES_match{imatch}(1) ;
    NODES_i = MESH.NODES_FACES{iface} ; 
    f{iface} = small2large(NODES_i,ndim) ; 
    COOR_i = MESH.PROPERTIES_FACES{iface}.COORrelA_global ;  
    
    
    jface = FACES_match{imatch}(2) ;
    COOR_j = MESH.PROPERTIES_FACES{jface}.COORrelA_global ; 
    NODES_j = MESH.NODES_FACES{jface} ; 
    f{jface} = small2large(NODES_j,ndim) ; 
    Idx_j = knnsearch(COOR_j,COOR_i) ;
    
    COOR_j_i = COOR_j(Idx_j,:) ;  
    
    
    errorCOOR = norm(COOR_j_i-COOR_j,'fro') ;
    disp('---------------------------------------------')
    disp(['Checking matching interfaces = ',num2str(iface),' and ',num2str(jface)]) ;
    disp(['ERROR_match = ',num2str(errorCOOR)])
    if errorCOOR > 1e-6
        error('No periodicity in this problem')
    else
       % NODES_j(Idx_j)
          Idx_DOFs_j_i = small2large(Idx_j,ndim) ;
         f{jface} = f{jface}(Idx_DOFs_j_i) ;
    end
    
    
     
     % Local numbering (with respect to b)
    [dummy1,F{iface},dummy2] =intersect(b,f{iface},'stable');
    [dummy1,F{jface},dummy2] =intersect(b,f{jface},'stable');

    
end


% $\Z = \Mbar \R_b$
% -----------------------------
% Construction of Mbar (Geometric mass matrix, interface boundaries)
%      f = unique(cell2mat(faceNODES(:))) ;
ndimLOC = 1;
wDIAGb = CompWeightDiag(GEOproperties.wSTb,ndimLOC)  ;
Mbound = GEOproperties.NstB_left'*(wDIAGb*GEOproperties.NstB_left) ;
Mff = Mbound(rnodLOC,rnodLOC) ;
% This is a matrix relating scalar nodal vectors. Extension to vectors of
% ndim entries per node
nnodeB = size(Mff,1) ;
ndim = 2 ;
Mbar = sparse(nnodeB*ndim,nnodeB*ndim) ;
for idim = 1:ndim
    Mbar(idim:ndim:end,idim:ndim:end)  = Mff ;
end

Rb = GEOproperties.RIGID_BODY_MODES(b,:) ; 
%PROP_GEO = Rb'*Mbar*Rb ; 

%%%  where $\Z = \Mbar \R_b$, and $\Z_{F_i} = \Z(\F_i,:)$.  

Z = Mbar*Rb ; 

%  \item Let $\alphaB \subset \{ 1,2 \ldots M \}$ be a subset of 3 entries of $\F_3$ so that 
% $\Z(\F_3(\alpha),:)$ is invertible, and let $\betaB$ denote the complementary set in $\F_3$. 
% With these decomposition at hand, \refeq{eq:4,dass} may be written as 

 [~,alphaB]=licols(Z(F{3},:)') ; %
 % Complementary set 
 betaB = setdiff(1:length(F{3}),alphaB) ; 
