function [MESHcluster] =AssignPropSlice3D_1D(MESH1D,MESH3D)
% Inputs: MESH1D --> Structure array containing data 1D skeleton
% MESH3D --> Structure array containint 3D data slice
% iacumCN --> Index of accumulated
% Outputs
% Global coordinates and connectivities (for printing purpose)
% Rotation matrices  (rotMAT)
% JAHO, 10-January-2018
% ----------------------

itype = MESH3D.INDEX ; % Index associated with the slice
IND_ELEM_TYPE = find(MESH1D.MaterialType == itype) ;  % Domains pertaining to type "itype"
CNintLOC = MESH1D.CN(IND_ELEM_TYPE,:) ; % 1D connectivities associated to this slice
nDOM = length(IND_ELEM_TYPE) ;   % Number of domains
SLICE = MESH3D.DATA3D;  % Mesh information Slice under study

% Angle of rotation around the x-axis
% The same for all except for providing domainwise information
% So far this option has not been implemented. Accordingly
MESH3D = DefaultField(MESH3D,'ROTATION_AROUND_LOCALaxis',0) ;
AngleRotationXaxis =  MESH3D.ROTATION_AROUND_LOCALaxis*ones(size(IND_ELEM_TYPE)) ;
%

ndim = size(SLICE.COOR,2) ; % Spatial dimensions, 2 or 3
coorINT = MESH1D.COOR ;
if size(coorINT,2) < ndim & itype == 1
    coorINT = [coorINT, zeros(size(coorINT,1),1)] ;
end
% Distance between centroids
DIST_CENTR = norm(SLICE.CENTRf2-SLICE.CENTRf1);
nelem = size(SLICE.CN,1) ; % Number of elements per slice
nnodeE = size(SLICE.CN,2) ; % Number of nodes per element
nnode = size(SLICE.COOR,1) ; % Number of nodes per slice
nnodeGLO = nDOM*nnode ;   % Total number of nodes
nelemGLO = nDOM*nelem ;  % Total number of elements
TypeElem = SLICE.TypeElement ; % Element type
COORglo = zeros(ndim,nnodeGLO ) ; % Global coordinates
CNglo = zeros(nelemGLO,nnodeE ) ; % Global coordinates
MaterialTypeglo = zeros(nelemGLO,1) ; % Global coordinates

rotMAT = zeros(ndim,ndim*size(CNintLOC,1)) ; % Global rotation matrix


ifinCN = 0 ;ifinCOOR  = 0 ; ifin =0 ;
% if itype == 1
%     iacumCN =   size(coorINT,1) ;
%
% end
for e = 1:nDOM
    nodoINI = CNintLOC(e,1) ;
    nodoFIN = CNintLOC(e,2) ;
    xINI = coorINT(nodoINI,:) ;
    xFIN = coorINT(nodoFIN,:) ;
    DIST = norm(xFIN-xINI) ;
    if abs(DIST- DIST_CENTR)>1e-2*DIST_CENTR
        error(['The width of this elements is ',num2str(DIST_CENTR),' m' ])
    end
    % Rotation matrix
    r1 = (xFIN-xINI)'/DIST ;  % Local x vector expressed in the global axis
    % Vector perpendicular to the plane defined by r1 and the global x axis
    %
    if norm(r1-[1 0 0]') >0
        nr = cross(r1,[1 0 0]') ;
        % Second vector  is perpendicular to r1 and nr
        r2 = cross(r1,nr) ;
        r2 = r2/norm(r2);
        % Third vector
        r3 = cross(r1,r2) ;
        rotMATloc = [r1,r2,r3] ;
    else
        rotMATloc = eye(3) ;
    end
    
    
    
    %% Rotation around local x-axix
    if AngleRotationXaxis(e) ~=0
        ANG  = AngleRotationXaxis(e) ;
        rotX = [1    0       0
            0  cosd(ANG) -sind(ANG)
            0  sind(ANG) cosd(ANG)      ] ;
        rotMATloc = rotX*rotMATloc ;
    end
    
    %%%
    
    %%%
    iini =   ifin +1 ;
    ifin = iini + ndim -1;
    rotMAT(:,iini:ifin) = rotMATloc;
    %%%
    % Global coordinate matrix
    % -------------------------
    % Coordinates relative to centroid face 1
    COORrel  = zeros(size(SLICE.COOR')) ;
    for idim=1:ndim
        COORrel(idim,:) = SLICE.COOR(:,idim) - SLICE.CENTRf1(idim) ;
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
    CNglo(iniCN:ifinCN,:) = SLICE.CN   ;
    
    MaterialTypeglo(iniCN:ifinCN) = SLICE.MaterialType;
    iacumCN  = iacumCN + nnode ;
end

% rotMATglo{itype} = rotMAT ;
MESHcluster.CN = CNglo ; 
MESHcluster.COOR = COORglo ; 
MESHcluster.rotMAT = rotMAT ; 
MESHcluster.MaterialType = MaterialTypeglo ; 


