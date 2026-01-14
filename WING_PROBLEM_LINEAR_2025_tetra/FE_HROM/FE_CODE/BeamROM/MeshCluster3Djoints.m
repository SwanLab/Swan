function [MESHcluster] =MeshCluster3Djoints(MESH1D,MESH3D,ientity,itypeJOINT,DATA)
% Clustering of mesh data (joints )
% Inputs: MESH1D --> Structure array containing data 1D skeleton
% MESH3D --> Structure array containint 3D data slice
% Outputs
% Global coordinates and connectivities
% Rotation matrices  (rotMAT)
% JAHO, 9-FEbruary-2018
% ----------------------
if nargin == 0
    load('tmp.mat')
end

JOINT = MESH3D.JOINTS(itypeJOINT).DATA3D;  % Mesh information JOINT under study
iacumCN = 0 ;
ndim = size(JOINT.COOR,2) ; % Spatial dimensions, 2 or 3
nelem = size(JOINT.CN,1) ; % Number of elements per JOINT

JOINT = DefaultField(JOINT,'CNb',[]) ;
if isempty(JOINT.CNb)
    JOINT.CNb = JOINT.CONNECTb ;
end
nelemB = size(JOINT.CNb,1) ; % Number of bound. elements per JOINT

nnodeE = size(JOINT.CN,2) ; % Number of nodes per element
nnodeEb = size(JOINT.CNb,2) ; % Number of nodes per element

nnode = size(JOINT.COOR,1) ; % Number of nodes per JOINT

%IND_BEAMS= find(MESH1D.Elements2Print(:,2) == ientity) ;  % Domains pertaining to type "ientity"

%IND_JOINTS = find(MESH1D.Elements2Print(IND_BEAMS,3)==itypeJOINT) ;
%IND_JOINTS =IND_BEAMS(IND_JOINTS) ;

%IND_ELEM_TYPE = MESH1D.Elements2Print(IND_JOINTS,1) ;

%nDOM = length(IND_ELEM_TYPE) ;   % Number of domains
%coorINT = MESH1D.COOR ;

% END_NODES  joints of this type 
END_NODES = MESH1D.INFOLINES.END_NODES{ientity} ; 
END_ELEMENTS = MESH1D.INFOLINES.END_ELEMENTS{ientity} ; 

nDOM = length(END_NODES) ; % Number of domains

nnodeGLO = nDOM*nnode ;   % Total number of nodes
nelemGLO = nDOM*nelem ;  % Total number of elements
COORglo = zeros(ndim,nnodeGLO ) ; % Global coordinates
CNglo = zeros(nelemGLO,nnodeE ) ; % Global con.
if iscell(JOINT.CNb)
    CNbGLO = [] ; %cell(size(SLICE.CNb)) ;
else
   nelemBglo = nDOM*nelemB ;  % Total number of bound. elements

   CNbGLO = zeros(nelemBglo,nnodeEb ) ; % Global con. bound.

end
MaterialTypeglo = zeros(nelemGLO,1) ; %
ifinCN = 0 ;ifinCOOR  = 0 ;  ifinCNb = 0 ;


for ijoint = 1:nDOM   % Loop over joints 
    
    nodoINI = END_NODES{ijoint}(1) ; % Centroid face 1     
    elemINI = END_ELEMENTS{ijoint}(1) ;
    xINI = MESH1D.COOR(nodoINI,:) ; % Coordinates
    
    % Rotation matrix
    ifinROT = ndim*elemINI ;
    iiniROT = ifinROT-ndim+1 ;
    rotMATloc  = MESH1D.ROTATIONS(:,iiniROT:ifinROT)  ; 
    % Global coordinate matrix
    % -------------------------
    % Coordinates relative to centroid face 1
    COORrel  = zeros(size(JOINT.COOR')) ;
    
    if isfield(JOINT,'CENTRf1')
    CF1 = JOINT.CENTRf1 ; 
    else
        CF1 = JOINT.CENTROIDS(:,1); 
    end
    
    for idim=1:ndim
        COORrel(idim,:) = JOINT.COOR(:,idim) - CF1(idim) ;
    end
  
    % Rotated coordinates
    COORrel = rotMATloc*COORrel ;
    % Translation
    for idim = 1:ndim
        COORrel(idim,:) = COORrel(idim,:) + xINI(idim) ;
    end
    % ------------
    iiniCO0R = ifinCOOR + 1;
    ifinCOOR = iiniCO0R + nnode -1;
    COORglo(:,iiniCO0R:ifinCOOR) = COORrel ;
    
    % Global connectivities
    iniCN = ifinCN + 1;
    ifinCN = iniCN + nelem -1 ;
    CNglo(iniCN:ifinCN,:) = JOINT.CN +  iacumCN ;
    % Bnd. connectivities
   % iniCNb = ifinCNb + 1;
   % ifinCNb = iniCNb + nelemB -1 ;
   % CNbGLO(iniCNb:ifinCNb,:) = JOINT.CNb +  iacumCN ;
    
    MaterialTypeglo(iniCN:ifinCN) = JOINT.MaterialType;
    iacumCN  = iacumCN + nnode ;
end

MESHcluster.CN = CNglo; ;
MESHcluster.CNb = CNbGLO; ;

MESHcluster.COOR = COORglo' ;
MESHcluster.MaterialType = MaterialTypeglo ;


