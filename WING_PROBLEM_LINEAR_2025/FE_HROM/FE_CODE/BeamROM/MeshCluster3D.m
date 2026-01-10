function [MESHcluster] =MeshCluster3D(MESH1D,MESH3D,iENTITY,itypeSLICE,DATAIN)
% Clustering of mesh data (slices )
% Inputs: MESH1D --> Structure array containing data 1D skeleton
% MESH3D --> Structure array containint 3D data slice
% Outputs
% Global coordinates and connectivities
% Rotation matrices  (rotMAT)
% JAHO, 10-January-2018
% ----------------------
if nargin == 0
    load('tmp2.mat')
    
end

SLICE = MESH3D.SLICES(itypeSLICE).DATA3D;  % Mesh information Slice under study
iacumCN = 0 ;
ndim = size(SLICE.COOR,2) ; % Spatial dimensions, 2 or 3

DATAIN = DefaultField(DATAIN,'setReducedElements',[]) ; 

if    ~isempty(DATAIN.setReducedElements)
    SLICE.CN = SLICE.CN(DATAIN.setReducedElements,:) ; % Only reduced set of elements 
    SLICE.MaterialType = SLICE.MaterialType(DATAIN.setReducedElements,:) ;
end

nelem = size(SLICE.CN,1) ; % Number of elements per slice


SLICE = DefaultField(SLICE,'CNb',[]) ;
if isempty(SLICE.CNb)
    SLICE.CNb = SLICE.CONNECTb ;
end
nelemB = size(SLICE.CNb,1) ; % Number of bound. elements per slice

nnodeE = size(SLICE.CN,2) ; % Number of nodes per element
nnodeEb = size(SLICE.CNb,2) ; % Number of nodes per element

nnode = size(SLICE.COOR,1) ; % Number of nodes per slice

IND_BEAMS = find(MESH1D.MaterialType == iENTITY) ; % Domains pertaining to beam type "iENTITY". Global numbering

%MESH1D = DefaultField(MESH1D,'TypeOfSlices',ones(size(MESH1D.MaterialType))) ;
%IND_SLICES = find(MESH1D.TypeOfSlices(IND_BEAMS) == itypeSLICE) ; % Domains pertaining to slice type "itypeSLICE"
% Intersection between these two sets
IND_SLICES_INTER = IND_BEAMS;
% Indices to print
DATAIN = DefaultField(DATAIN,'SELECTED_DOMAINS',[]) ;

if isempty(DATAIN.SELECTED_DOMAINS)
    MESH1D = DefaultField(MESH1D,'Elements2Print',IND_SLICES_INTER) ;
    IND_ELEM_TYPE = intersect(IND_SLICES_INTER,MESH1D.Elements2Print) ;
else
    %  IND_ELEM_TYPE =  (DATAIN.SELECTED_DOMAINS{iENTITY}) ;
    IND_ELEM_TYPE = IND_BEAMS(DATAIN.SELECTED_DOMAINS{iENTITY}) ;
end

%  LOCAL_INDEX
IND_LOCAL_3D = zeros(size(IND_ELEM_TYPE)) ;
for  ilocal = 1:length(IND_ELEM_TYPE)
    IND_LOCAL_3D(ilocal) =    find(IND_ELEM_TYPE(ilocal) == IND_SLICES_INTER) ;
end

%%  Lateral surfaces
% --------------------
PROPS  =MESH1D.PROP(iENTITY)  ;
PROPS= DefaultField(PROPS,'LATERAL_SURFACES',[] ) ; %.PRINT_ALL
PROPS.LATERAL_SURFACES = ...
    DefaultField(PROPS.LATERAL_SURFACES,'PRINT_ALL',[] );

%

if PROPS.LATERAL_SURFACES.PRINT_ALL == 1
    
%     if ~isempty(DATAIN.setReducedElements)
%         error('Non compatible options (lateral surfaces and integration elements)')
%     end
    IND_ELEM_TYPE_INNER = IND_ELEM_TYPE ;  % Interior elements, 3D
    IND_ELEM_TYPE = IND_SLICES_INTER ;
    WHAT_TO_PRINT = zeros(size(IND_ELEM_TYPE)) ; % If = 1, then 3D
    %  WHAT_TO_PRINT([1,end]) = 1;
    WHAT_TO_PRINT(IND_LOCAL_3D) = 1;
    
    nDOM_3D = length(find(WHAT_TO_PRINT>0)) ;
    nDOM_lat = length(find(WHAT_TO_PRINT==0)) ;
    
    
    % Connectivities lateral surfaces
    if DATAIN.POST_PROCESS_LATERAL_SURFACES.COARSE_MESH == 1  & ~isempty(SLICE.CONNECTb_coarse)
        CNb_lateral  =  cell2mat(SLICE.CONNECTb_coarse(3:end)') ;
    else
        CNb_lateral  =  cell2mat(SLICE.CONNECTb(3:end)') ;
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
coorINT = MESH1D.COOR ;
% Distance between centroids
nnodeGLO = nDOM_3D*nnode + nDOM_lat*nnodesLATERAL ;   % Total number of nodes
nelemGLO_3D = nDOM_3D*nelem ;  % Total number of elements (3D)
nelemGLOlat= nDOM_lat*nelemL ;  % Total number of elements (lateral)

nelemGLO = nelemGLO_3D + nelemGLOlat ;

COORglo = zeros(ndim,nnodeGLO ) ; % Global coordinates
CNglo = zeros(nelemGLO_3D,nnodeE ) ; % Global conn, 3D
if   nelemGLOlat ~=0
    CNgloLAT = zeros(nelemGLOlat,nnodeBND ) ; % Global conn, lateral surfaces
else
    CNgloLAT = [] ;
end


if iscell(SLICE.CNb)
    CNbGLO = [] ; %cell(size(SLICE.CNb)) ;
else
    nelemBglo = nDOM*nelemB ;  % Total number of bound. elements
    CNbGLO = zeros(nelemBglo,nnodeEb ) ; % Global con. bound.
end


MaterialTypeglo = zeros(nelemGLO_3D,1) ; %
MaterialTypeglo_lat = zeros(nelemGLOlat,1) ; %

ifinCN = 0 ;   ifinCN_lat = 0 ;   ; ifinCOOR  = 0 ;  ifinCNb = 0 ;

PROPS = MESH1D.PROP(iENTITY) ;
PROPS = DefaultField(PROPS,'xLOCAL',1) ;

DATAIN = DefaultField(DATAIN,'xSIGN',[]) ; 

if ~isempty(DATAIN.xSIGN)
    SIGNxLOCAL = DATAIN.xSIGN ; 
    % Amendment, 12-Dec-2018... To be revised, not general
else
    if isempty(PROPS.xLOCAL)
        SIGNxLOCAL = 1;
    else
        SIGNxLOCAL =PROPS.xLOCAL;
    end
end

MESH1D  = DefaultField( MESH1D,'TRANSFM',[]) ;

DATAIN = DefaultField(DATAIN,'SEPARATION_BETWEEN_SLICES_PLOT',zeros(1,ndim)) ;

for e = 1:nDOM
    eGLO = IND_ELEM_TYPE(e) ;
    if SIGNxLOCAL == 1
        nodoINI = MESH1D.CN(eGLO,1) ;
        nodoFIN = MESH1D.CN(eGLO,2) ;
    else
        nodoINI = MESH1D.CN(eGLO,2) ;
        nodoFIN = MESH1D.CN(eGLO,1) ;
    end
    
    xINI = MESH1D.COOR(nodoINI,:) ;
    xFIN = MESH1D.COOR(nodoFIN,:) ;
  %  DIST = norm(xFIN-xINI) ;
    %DIST_SEPAR = DIST*DATAIN.SEPARATION_BETWEEN_SLICES_PLOT ;
    xINI_plot = xINI ; %+ DIST_SEPAR*(e-1) ;
    
    %     if abs(DIST- DIST_CENTR)>1e-2*DIST_CENTR
    %         error(['The width of this elements is ',num2str(DIST_CENTR),' m' ])
    %     end
    % Rotation matrix
    ifinROT = ndim*eGLO ;
    iiniROT = ifinROT-ndim+1 ;
    rotMATloc  = MESH1D.ROTATIONS(:,iiniROT:ifinROT)  ;
    % TRANSFORMATION MATRICES
    if ~isempty(MESH1D.TRANSFM)
        Aloc = MESH1D.TRANSFM.A(:,iiniROT:ifinROT) ;
        Dloc = MESH1D.TRANSFM.D(:,iiniROT:ifinROT) ;
        a0loc = MESH1D.TRANSFM.a0(:,eGLO) ;
        
    end
    
    
    % Global coordinate matrix
    % -------------------------
    % Coordinates relative to centroid face 1
    COORrel  = zeros(size(SLICE.COOR')) ;
    for idim=1:ndim
        COORrel(idim,:) = SLICE.COOR(:,idim) - SLICE.CENTRf1(idim) ;
    end
    
    if ~isempty(MESH1D.TRANSFM)
        % Transformation (scaling, varying cross-section)
        % x = a0 + A*X +X_1*D*X
        COORtransf = zeros(size(COORrel)) ;
        for idim = 1:ndim
            COORtransf(idim,:) = a0loc(idim) + Aloc(idim,:)*COORrel + COORrel(1,:).*(Dloc(idim,:)*COORrel) ;
        end
        
    else
        COORtransf = COORrel ;
    end
    
    % Rotated coordinates
    COORrel = rotMATloc*COORtransf ;
    % Translation
    for idim = 1:ndim
        COORrel(idim,:) = COORrel(idim,:) + xINI_plot(idim) ;
    end
    % ------------
    
    
    % Global connectivities
    if WHAT_TO_PRINT(e) == 1
        
        iiniCO0R = ifinCOOR + 1;
        ifinCOOR = iiniCO0R + nnode -1;
        COORglo(:,iiniCO0R:ifinCOOR) = COORrel ;
        
        iniCN = ifinCN + 1;  % 3D slices
        ifinCN = iniCN + nelem -1 ;
        CNglo(iniCN:ifinCN,:) = SLICE.CN +  iacumCN ;
        MaterialTypeglo(iniCN:ifinCN) = SLICE.MaterialType;
        iacumCN  = iacumCN + nnode ;
    else
        
        %         NODESlateral  = unique(CNb_lateral) ;
        %     nnodesLATERAL = length(NODESlateral) ;
        iiniCO0R = ifinCOOR + 1;
        ifinCOOR = iiniCO0R + nnodesLATERAL -1;
        COORglo(:,iiniCO0R:ifinCOOR) = COORrel(:,NODESlateral) ;
        
        % Lateral surfaces
        iniCN_lat = ifinCN_lat + 1;  % 3D slices
        ifinCN_lat = iniCN_lat + nelemL -1 ;
        CNgloLAT(iniCN_lat:ifinCN_lat,:) = CNb_lateral +  iacumCN ;
        MaterialTypeglo_lat(iniCN_lat:ifinCN_lat) = (max(SLICE.MaterialType)+1)*ones(size(iniCN_lat:ifinCN_lat));
        iacumCN  = iacumCN + nnodesLATERAL ;
        
    end
    
    
    %     % Bnd. connectivities
    %     iniCNb = ifinCNb + 1;
    %     ifinCNb = iniCNb + nelemB -1 ;
    %     if iscell(SLICE.CNb)
    %
    %     else
    %         CNbGLO(iniCNb:ifinCNb,:) = SLICE.CNb +  iacumCN ;
    %     end
    
    
end

MESHcluster.CN = CNglo; ;
MESHcluster.CNb = CNbGLO; ;

MESHcluster.CNlateral = CNgloLAT; ;
MESHcluster.CNb = CNbGLO; ;

MESHcluster.COOR = COORglo' ;
MESHcluster.MaterialType = MaterialTypeglo ;
MESHcluster.MaterialTypeglo_lat = MaterialTypeglo_lat ;


