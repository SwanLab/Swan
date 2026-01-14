function [MESHcluster] =MeshCluster3D_RVE(MESH2D,MESH3D,iENTITY,itypeSLICE,DATAIN)
% Clustering of mesh data (slices )
% Inputs: MESH2D --> Structure array containing data 2D skeleton
% MESH3D --> Structure array containint 3D data slice
% Outputs
% Global coordinates and connectivities
% JAHO, 20-July-2018
% ----------------------
if nargin == 0
    load('tmp1.mat')
    
end

RVE = MESH3D.RVES(itypeSLICE).DATA3D;  % Mesh information Slice under study
iacumCN = 0 ;
ndim = size(RVE.COOR,2) ; % Spatial dimensions, 2 or 3
nelem = size(RVE.CN,1) ; % Number of elements per slice
RVE = DefaultField(RVE,'CNb',[]) ;
if isempty(RVE.CNb)
    RVE.CNb = RVE.CONNECTb ;
end
nelemB = size(RVE.CNb,1) ; % Number of bound. elements per RVE

nnodeE = size(RVE.CN,2) ; % Number of nodes per element
nnodeEb = size(RVE.CNb,2) ; % Number of nodes per element

nnode = size(RVE.COOR,1) ; % Number of nodes per RVE
nfaces_CONNECT = 4 ; 

IND_RVES = find(MESH2D.MaterialType == iENTITY) ; % Domains pertaining to beam type "iENTITY". Global numbering

% Intersection between these two sets
IND_RVES_INTER = IND_RVES;
% Indices to print
DATAIN = DefaultField(DATAIN,'SELECTED_DOMAINS',[]) ;

if isempty(DATAIN.SELECTED_DOMAINS)
    MESH2D = DefaultField(MESH2D,'Elements2Print',IND_RVES_INTER) ;
    IND_ELEM_TYPE = intersect(IND_RVES_INTER,MESH2D.Elements2Print) ;
else
    %  IND_ELEM_TYPE =  (DATAIN.SELECTED_DOMAINS{iENTITY}) ;
    IND_ELEM_TYPE = IND_RVES(DATAIN.SELECTED_DOMAINS{iENTITY}) ;
end

%  LOCAL_INDES
IND_LOCAL_3D = zeros(size(IND_ELEM_TYPE)) ;
for  ilocal = 1:length(IND_ELEM_TYPE)
    IND_LOCAL_3D(ilocal) =    find(IND_ELEM_TYPE(ilocal) == IND_RVES_INTER) ;
end

%%  Lateral surfaces
% --------------------
PROPS  =MESH2D.PROP(iENTITY)  ;
PROPS= DefaultField(PROPS,'LATERAL_SURFACES',[] ) ; %.PRINT_ALL
PROPS.LATERAL_SURFACES = ...
    DefaultField(PROPS.LATERAL_SURFACES,'PRINT_ALL',[] );

%

if PROPS.LATERAL_SURFACES.PRINT_ALL == 1
    IND_ELEM_TYPE_INNER = IND_ELEM_TYPE ;  % Interior elements, 3D
    IND_ELEM_TYPE = IND_RVES_INTER ;
    WHAT_TO_PRINT = zeros(size(IND_ELEM_TYPE)) ; % If = 1, then 3D
    %  WHAT_TO_PRINT([1,end]) = 1;
    WHAT_TO_PRINT(IND_LOCAL_3D) = 1;
    
    nDOM_3D = length(find(WHAT_TO_PRINT>0)) ;
    nDOM_lat = length(find(WHAT_TO_PRINT==0)) ;
    
    
    % Connectivities lateral surfaces
    if DATAIN.POST_PROCESS_LATERAL_SURFACES.COARSE_MESH == 1  & ~isempty(RVE.CONNECTb_coarse)
        CNb_lateral  =  cell2mat(RVE.CONNECTb_coarse(nfaces_CONNECT+1:end)') ;
    else
        CNb_lateral  =  cell2mat(RVE.CONNECTb(nfaces_CONNECT+1:end)') ;
    end
    
    NODESlateral  = unique(CNb_lateral) ;
    nnodesLATERAL = length(NODESlateral) ;
    CNb_lateral = RenumberConnectivities(CNb_lateral,1:nnodesLATERAL);
    
    
    nelemL= size(CNb_lateral,1) ;
    nnodeBND = size(CNb_lateral,2) ;
    
else
    CNb_lateral = [] ;
    nelemGLOlat = 0 ;
    nelemL = 0 ;
    nDOM_3D = length(IND_ELEM_TYPE) ;
    nDOM_lat = 0 ;
    nnodesLATERAL = 0 ;
    WHAT_TO_PRINT = ones(size(IND_ELEM_TYPE) ) ;
end


nDOM = nDOM_3D+nDOM_lat ;   % Number of domains
% Distance between centroids
nnodeGLO = nDOM_3D*nnode + nDOM_lat*nnodesLATERAL ;   % Total number of nodes
nelemGLO_3D = nDOM_3D*nelem ;  % Total number of elements (3D)
nelemGLOlat= nDOM_lat*nelemL ;  % Total number of elements (lateral)

 
COORglo = zeros(ndim,nnodeGLO ) ; % Global coordinates
CNglo = zeros(nelemGLO_3D,nnodeE ) ; % Global conn, 3D
if   nelemGLOlat ~=0
    CNgloLAT = zeros(nelemGLOlat,nnodeBND ) ; % Global conn, lateral surfaces
else
    CNgloLAT = [] ;
end


if iscell(RVE.CNb)
    CNbGLO = [] ; %cell(size(RVE.CNb)) ;
else
    nelemBglo = nDOM*nelemB ;  % Total number of bound. elements
    CNbGLO = zeros(nelemBglo,nnodeEb ) ; % Global con. bound.
end


MaterialTypeglo = zeros(nelemGLO_3D,1) ; %
MaterialTypeglo_lat = zeros(nelemGLOlat,1) ; %

ifinCN = 0 ;   ifinCN_lat = 0 ;   ; ifinCOOR  = 0 ;  ifinCNb = 0 ;

 


 
for e = 1:nDOM
    eGLO = IND_ELEM_TYPE(e) ;
    
    % Relative coordinates of the RVE  with respect to the CORNER_1_CENTROID
    % face 1
    COORrel  = zeros(size(RVE.COOR')) ;
    for idim=1:ndim
        COORrel(idim,:) = RVE.COOR(:,idim) - RVE.CORNER_1_CENTROID(idim) ;
    end
    
    % CORNER_1_CENTROID corresponds to the first node 
    nodoINI = MESH2D.CN(eGLO,1) ; 
    cOORloc = MESH2D.COOR(nodoINI,:) ; 
  %  [~,IIII] = min(cOORloc(:,1)) ;  
     
    xINI = cOORloc  ; 
    % Translation
    for idim = 1:ndim
        COORrel(idim,:) = COORrel(idim,:) + xINI(idim) ;
    end
    % ------------    
    
    % Global connectivities
    if WHAT_TO_PRINT(e) == 1
        
        iiniCO0R = ifinCOOR + 1;
        ifinCOOR = iiniCO0R + nnode -1;
        COORglo(:,iiniCO0R:ifinCOOR) = COORrel ;
        
        iniCN = ifinCN + 1;  % 3D RVEs
        ifinCN = iniCN + nelem -1 ;
        CNglo(iniCN:ifinCN,:) = RVE.CN +  iacumCN ;
        MaterialTypeglo(iniCN:ifinCN) = RVE.MaterialType;
        iacumCN  = iacumCN + nnode ;
    else
        
        %         NODESlateral  = unique(CNb_lateral) ;
        %     nnodesLATERAL = length(NODESlateral) ;
        iiniCO0R = ifinCOOR + 1;
        ifinCOOR = iiniCO0R + nnodesLATERAL -1;
        COORglo(:,iiniCO0R:ifinCOOR) = COORrel(:,NODESlateral) ;
        
        % Lateral surfaces
        iniCN_lat = ifinCN_lat + 1;  % 3D RVEs
        ifinCN_lat = iniCN_lat + nelemL -1 ;
        CNgloLAT(iniCN_lat:ifinCN_lat,:) = CNb_lateral +  iacumCN ;
        MaterialTypeglo_lat(iniCN_lat:ifinCN_lat) = (max(RVE.MaterialType)+1)*ones(size(iniCN_lat:ifinCN_lat));
        iacumCN  = iacumCN + nnodesLATERAL ;
        
    end
    
    
    %     % Bnd. connectivities
    %     iniCNb = ifinCNb + 1;
    %     ifinCNb = iniCNb + nelemB -1 ;
    %     if iscell(RVE.CNb)
    %
    %     else
    %         CNbGLO(iniCNb:ifinCNb,:) = RVE.CNb +  iacumCN ;
    %     end
    
    
end

MESHcluster.CN = CNglo; ;
MESHcluster.CNb = CNbGLO; ;

MESHcluster.CNlateral = CNgloLAT; ;
MESHcluster.CNb = CNbGLO; ;

MESHcluster.COOR = COORglo' ;
MESHcluster.MaterialType = MaterialTypeglo ;
MESHcluster.MaterialTypeglo_lat = MaterialTypeglo_lat ;


