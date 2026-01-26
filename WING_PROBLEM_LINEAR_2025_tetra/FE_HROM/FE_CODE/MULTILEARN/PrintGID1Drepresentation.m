function PrintGID1Drepresentation(DATAINM,dFACES,DATAONLINE,pINT,alphaBC,BasisINT,COORref,NODESfaces,...
    BasisRdef,BasisRrb,rDEF,rRB,MAXstressVONMISES)

if nargin == 0
    load('tmp.mat')
end

if DATAINM.MODES_INTERFACE_RIGID_BODY ~=1
    error('Option not compatible')
end

[DATA,COOR,CN,d] = Construct1D_mesh_FORCES(DATAINM,dFACES,DATAONLINE,pINT,alphaBC,BasisINT,COORref,NODESfaces,...
    BasisRdef,BasisRrb,rDEF,rRB); 

strainGLO =[] ;
stressGLO = [] ;
TypeElement = 'Linear' ;
React = [] ;
posgp = [] ;
NAME_INPUT_DATA = 'Beam1D' ;
NameFileMesh = '' ;
MaterialType = [];
DATA.MAXstressVONMISES = MAXstressVONMISES ;

GidPostProcess(COOR,CN,TypeElement,d,strainGLO, stressGLO, ...
    React,NAME_INPUT_DATA,posgp,NameFileMesh,MaterialType,DATA);
