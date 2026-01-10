function PrintGID1D_3Drepresentation(DATAINM,dFACES,DATAONLINE,pINT,alphaBC,BasisINT,COORref,NODESfaces,...
    BasisRdef,BasisRrb,rDEF,rRB,MAXstressVONMISES,...
    COOR_3d,CN_3d,TypeElement_3d,d_3d, strainGLO_3d,stressGLO_3d,  REACT_3d,NAME_BASE_GIDfiles, ...
    posgp_3d,NameFileMeshs_3d,MaterialType,DATA,RESIDUAL)

if nargin == 0
    load('tmp.mat')
end

DATA = DefaultField(DATA,'PostProcessWithNOSLICES',1) ; 
DATA = DefaultField(DATA,'MATERIAL',[]) ; 

if DATA.PostProcessWithNOSLICES == 1 & ~isempty(DATA.MATERIAL)
    nmat = length(DATA.MATERIAL.PLY) ; 
    ndom = DATAONLINE.NdomX ; 
    TOTnmat = nmat*ndom ; 
    MAT_TypeORIG = (1:TOTnmat) ;
    MAT_TypeNOSLICES = repmat(((1:nmat)'+1),ndom,1) ; 
    NewMaterialType = zeros(size(MaterialType)); 
   
    for imat = 1:length(MAT_TypeORIG)
        INDLOC =find(MaterialType==imat) ; 
        NewMaterialType(INDLOC) = MAT_TypeNOSLICES(imat) ; 
    end
    
   MaterialType =  NewMaterialType ; 
end


% if DATAINM.MODES_INTERFACE_RIGID_BODY ~=1
%     error('Option not compatible')
% end

[DATA,COOR_1d,CN_1d,d_1d] = Construct1D_mesh_FORCES(DATAINM,dFACES,DATAONLINE,pINT,alphaBC,BasisINT,COORref,NODESfaces,...
    BasisRdef,BasisRrb,rDEF,rRB,DATA); 

COOR = [COOR_1d; COOR_3d]; % COORDINATES
NODES = {(1:size(COOR_1d,1))', (1:size(COOR_3d,1))'+size(COOR_1d,1)}; 
displacement = [d_1d; d_3d] ; % DISPLACEMENTS
REACT = [] ; % 
 
% CONNECTIVITIES 
CN = {CN_1d,CN_3d+size(COOR_1d,1)} ; 
DATA.NAMEMESH = {'MESH_1D','MESH_3D'} ; 
MaterialType = {ones(size(CN_1d,1),1),MaterialType+1} ; 
TypeElement = {'Linear',TypeElement_3d };

if DATAINM.PRINT_AVERAGE_STRESSES_ON_ELEMENTS ==1
    posgp_3d = [] ; 
end
posgp = {[],posgp_3d} ; 

strainGLO =[] ;

React = [] ;NAME_INPUT_DATA = 'Beam1D_3D' ;
NameFileMesh = '' ;
DATA.MAXstressVONMISES = MAXstressVONMISES ;

GidPostProcess_multi(COOR,CN,TypeElement,displacement,strainGLO, stressGLO_3d, ...
    React,NAME_INPUT_DATA,posgp,NameFileMesh,MaterialType,DATA,NODES,RESIDUAL,DATAINM);
