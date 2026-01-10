function [MESH1D,xSIGN] = Geometry1Dstructure(MESH1Dinput,MESH3D,DATAIN)


% -----------------------------------------
% Reading coordinate structure (1D skeleton, linear)
% -----------------------------------------
[MESH1D]= ReadMeshFileStr(MESH1Dinput.NAME,'READ_MATERIAL_COLUMN',1)  ;
MESH1D.NAME = MESH1Dinput.NAME ;
MESH1D.PROP = MESH1Dinput.PROP ;
if size(MESH1D.COOR,2) ==2 && size(MESH3D.SLICES(1).DATA3D.COOR,2)==3
    MESH1D.COOR = [MESH1D.COOR,zeros(size(MESH1D.COOR,1),1)] ;
end

% READING .dat files (faces and lines).  
[MESH1D.NODES_POINTS,MESH1D.NODES_LINES] = NodeLinesPointsGID(MESH1Dinput.NAME) ; 
%%%%%%
 
% Mixed beams. Determine which type of slice is associated to each 1D
% element.  MESH1D.MaterialType specifies the type of beam, but beams may contain several types of
% slices. Basically, it assigns to each element their corresponding  3D
% index (either beam or joint)
MESH1D.TypeOfSlices = WhichTypeOfSlices(MESH1D) ;   % All = 1 , Eliminated. 29-May-2018- Homogeneous beams

% ------------------------------------------------------------
% SEPARATION OF ENTITIES DEFINED IN THE GID'S INPUT PROJECT - Disjoints
% structures
% -----------------------------------------------------------
[MESH1D.INFOLINES ]= SeparateBeamsAndJoints(MESH1D) ;
% --------------------------------------------------------
% DEFINITION LOCAL AXES 1d SKELETON. PRINTING LOCAL AXES (for beams !!!!)
% ---------------------------------------------------------
[MESH1D.ROTATIONS,MESH1D.ORDER_CONNECTIVITIES,xSIGN] = LocalAxisDefinitionBeams(MESH1D,MESH3D,DATAIN) ;

 

% --------------------------------------------------------
% DEFINITION TRANSFORMATION MATRICES (DILATATION and VARIATIONS OF CROSS-SECTIONAL AREA), for beams !
% ---------------------------------------------------------
MESH1D.TRANSFM = [] ;
DATAIN = DefaultField(DATAIN,'INTRODUCE_DILATAION_DOMAINS',0) ; 
if DATAIN.INTRODUCE_DILATAION_DOMAINS == 1
[MESH1D.TRANSFM.A,MESH1D.TRANSFM.D,MESH1D.TRANSFM.a0] = CrossSectionsVariations_dilat(MESH1D,MESH3D,DATAIN) ;
end

% Deciding which beam domains have to be printed
% --------------------------------------------
MESH1D.Elements2Print = WhichBeamElementsToPrint(MESH1D) ;