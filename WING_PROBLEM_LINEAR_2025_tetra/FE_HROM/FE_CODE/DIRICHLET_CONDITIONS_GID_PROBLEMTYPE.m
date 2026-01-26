% Prescribing general Dirichlet Boundary Conditions
%  JAHO.
%20-June-2019  , it is only valid for GLOBAL coordinates
% -----------------------------------------------------------------------------

% FIRST STEP ---> Lines and surfaces created with GID
% ---------------------------------------------------
% READING .dat files (faces and lines). Generated with GID's problemtype PROBLEMTYPES_GID/SLICE_JOINTS.gid
[NODES_FACES,NODES_LINES] = NodesFacesLinesGID(NameFileMesh) ;


if size(COOR,2)==2
    NODES_ENTITIES = NODES_LINES ;
    LABEL = 'LINE' ;
    nnn = min(length(NODES_ENTITIES),length(INPUTS_LOC.DISP.LINE)) ;
    
else
    NODES_ENTITIES = NODES_FACES ;
    LABEL = 'FACE' ;
    nnn = min(length(NODES_ENTITIES),length(INPUTS_LOC.DISP.FACE)) ;
    
end

 [Gb,dR,DOFr,DOFm] = DirichletConditionsNodesXYZ(INPUTS_LOC,LABEL,NODES_ENTITIES,nnn,...
     COOR,CONNECTb,TypeElementB,DATA) ; 