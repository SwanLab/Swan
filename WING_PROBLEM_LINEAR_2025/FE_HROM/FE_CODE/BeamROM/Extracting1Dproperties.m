function MESH1D = Extracting1Dproperties(MESH1Dinput,MESH3D,DATAIN)


% -----------------------------------------
% Reading coordinate structure (1D skeleton, linear)
% -----------------------------------------
[MESH1D]= ReadMeshFileStr(MESH1Dinput.NAME,'READ_MATERIAL_COLUMN',1)  ;
MESH1D.NAME = MESH1Dinput.NAME ;
MESH1D.PROP = MESH1Dinput.PROP ;
if size(MESH1D.COOR,2) ==2 % 3D
    MESH1D.COOR = [MESH1D.COOR,zeros(size(MESH1D.COOR,1),1)] ;
end
% ------------------------------------------------------------
% SEPARATION OF ENTITIES DEFINED IN THE GID'S INPUT PROJECT - Disjoints
% structures 
% -----------------------------------------------------------
MESH1D.INFOLINES = SeparateBeamsAndJoints(MESH1D) ; 

% --------------------------------------------------------
% DEFINITION LOCAL AXES 1d SKELETON. PRINTING LOCAL AXES (for beams !!!!) 
% ---------------------------------------------------------
[MESH1D.ROTATIONS, MESH1D.Elements2Print] = LocalAxisDefinitionBeams(MESH1D,MESH3D,DATAIN) ; 
% --------------------------------------------------------
% DEFINITION TRANSFORMATION MATRICES (DILATATION and VARIATIONS OF CROSS-SECTIONAL AREA), for beams !  
% ---------------------------------------------------------
MESH1D.TRANSFM = [] ; 
[MESH1D.TRANSFM.A,MESH1D.TRANSFM.D,MESH1D.TRANSFM.a0] = CrossSectionsVariations_dilat(MESH1D,MESH3D,DATAIN) ; 