function [dRVE,NODESfaces,reactDOM ]= ExtractDisplReactionLOC(INPUT_DATAFILE,DATA)
% 
% INCLUDE_MODES_RVE = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COORrve = {} ;
CNrve = {}  ;
dRVE = {}  ;
NODESrve = {}  ;
NODESfaces = {}  ;
RIGID_BODY_MOTIONglo = {} ;
CoordinatesChangeGLO = {} ;
reactDOM = {} ; % REactions of each domain

IndElements = {} ;
% -------------------------
% (EXTRACTION OF PROPERTIES FROM FE FILES)
% -------------------------
NAME_MOD4ES =[INPUT_DATAFILE,'.msh'];
outputFE = ['DATAWS/',INPUT_DATAFILE,'_WS.mat'] ;
eval(INPUT_DATAFILE);
load(outputFE,'MaterialType','COOR','CN','d','TypeElement','posgp','Nst','wSTs_RHS','fNOD','posgp_RHS',...
    'CNb','TypeElementB','Fpnt','Tnod','CONNECTb','stressGLO','wSTs','Bst') ;
% Number of materials/bodies
IND_MAT =unique(MaterialType) ;
% Number of original materiasl
nMAT = length(MATERIAL.PLY) ;
% Number of RVEs
nRVE = length(IND_MAT)/nMAT         ;
% Store the coordinates, connectivitites and displacememts of each cell in
% a cell array
ndim = size(COOR,2) ;
d = reshape(d,ndim,[])' ; % Displacement vector
% Tolerance for matching nodes %
TOL = ChooseTolerance(CN,COOR) ;
% --------------------------------------------
% Loop over number of DOMAINS
% ----------------------------------------------
for i  = 1:nRVE
 
    % ------------------------------------------------
    % IDENTIFICATION elements and nodes of domain "i"
    % ------------------------------------------------
    IndElementsLOC = [] ;
    for imat = 1:nMAT
        IMAT = (i-1)*nMAT + imat ;
        IndElementsLOC =[IndElementsLOC;  find(IND_MAT(IMAT) == MaterialType)] ;   % Elements of domain i
    end
    IndElements{i} = IndElementsLOC ;
    CNrve{i} =  CN(IndElements{i},:) ; % Connectivities of Domain i
    NODESrve{i} =unique( CNrve{i}(:)) ;  % Nodes of Domain i
    COORabs = COOR(NODESrve{i},:) ;  % Coordinates of DOmain i
    dRVE{i} = d( NODESrve{i},:) ;  % Displacements of domain i
    % -----------------------------------------------------------
    % Purging rotations  and translation from displacement vector (main ouput dRVE)
    %------------------------------------------------------------
    [NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,COORrve] ...
        = PurgingRotations(i,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
        COORrve);
    % ---------------------------------------------------------------
    % Computing self-equilibrated reactions for each domain
    % -------------------------------------------------------
    [reactDOMloc] = ReactionsRVE(i,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
        COORrve,Nst,wSTs_RHS,fNOD,NODESrve,CNrve,posgp,IndElements{i},posgp_RHS,...
        CNb,TypeElementB,Fpnt,Tnod,CONNECTb,COOR,stressGLO,wSTs,Bst,TypeElement);
    reactDOM{i} =  reactDOMloc  ;
end


NCELLS = length(dRVE) ;

%%% REFERENCE MESH.
refMESH = 1 ;
%% COORDINATES OF REFERENCE MESH
COORref = COORrve{refMESH} ;
COORrefORD = reshape(COORref',prod(size(COORref)),1) ;
% Now we re-order
COORnew = {} ;
ElementsREF = IndElements{refMESH} ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identification of similar nodes between domains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(COORrve)
    if i ~= refMESH
        IDX = knnsearch(COORrve{i},COORref) ;
        COORrve{i} =  COORrve{i}(IDX,:) ;
        CoordinatesChange=  CoordinatesChangeGLO{i}(IDX,:) ;
        CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChange)),1 );
        COORnew{i} =COORrefORD + CoordinatesChange ;
        
    %    RIGID_BODY_MOTIONglo{i} =  RIGID_BODY_MOTIONglo{i}(IDX,:) ;
    %    RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
        dRVE{i} =  dRVE{i}(IDX,:)  ;
        dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
        IDXdof = small2large(IDX,ndim) ;
        reactDOM{i} = reactDOM{i}(IDXdof) ;
        for j = 1:length(NODESfaces{i})
            NODESfacesALL = zeros(size(NODESrve{i})) ;
            NODESfacesALL(NODESfaces{i}{j}) = 1 ;
            NODESfacesALL = NODESfacesALL(IDX) ;
            NODESfaces{i}{j}  =find(NODESfacesALL == 1) ;
        end
        NODESrve{i} =  NODESrve{i}(IDX) ;
    else
        dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
        RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
        CoordinatesChange=  CoordinatesChangeGLO{i}  ;
        CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChangeGLO{i})),1 );
        COORnew{i} =COORrefORD + CoordinatesChange ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
