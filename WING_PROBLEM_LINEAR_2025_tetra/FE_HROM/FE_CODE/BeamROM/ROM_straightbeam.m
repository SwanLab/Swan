function ROM_straightbeam(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,NAME_LOAD_DATA,MESH1D,MESH3D,DATARUN)
if nargin == 0
    load('tmp2.mat')
end
% --------------------------------------
% Reduced-order model for a beam-like structure
% made by repeating 3D slices
% JAHO, 1-May/4-July-2018
% ----------------------------------------------
% Defining search paths
if ~exist('INPUTS_GENERAL','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/']);  end
if ~exist('SVD_dom','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/MULTILEARN/']);  end
if ~exist('GeometryStructure','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/BeamROM/']);  end
% ------------------
%
% 1D -DATA (Supporting, skeleton mesh) --> MESH1D
% ----------------------------------------
% TYPE OF STRUCTURAL ENTITY --> MESH1D.TYPE.INDEX_BEAM_JOINT,
% MESH1D.TYPE.ISBEAM, MESH1D.TYPE.INDEX3D
% ----------------------------------------------------------
NAME_LOAD_DATA_inp = NAME_LOAD_DATA;
NAME_LOAD_DATA = [FOLDER,filesep,NAME_LOAD_DATA]  ; % Input file
run(NAME_LOAD_DATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DISP = DefaultField(DISP,'LEFT_END',[] );

if ~isempty(DISP.LEFT_END)
    COOR1_x  = MESH1D.COOR(MESH1D.NODES_POINTS{1},1) ;
    COOR2_x  = MESH1D.COOR(MESH1D.NODES_POINTS{2},1) ;
    IND1 = 1;
    IND2 = 2;
    if COOR2_x < COOR1_x
        IND2 = 1;
        IND1 = 2 ;
    end
    MESH1D.LEFT_END_NODE = MESH1D.NODES_POINTS{IND1}; MESH1D.RIGHT_END_NODE = MESH1D.NODES_POINTS{IND2} ;
else
    % Inputs are introduced in a different way  (directly using points 1 and 2)
end






% Retrieving stiffness matrix and other reduced-order variables.
% -------------------------------------------------------------
TYPE_STRUCTURES = unique(MESH1D.MaterialType) ; % Number of different structural entities
NTYPE = length(TYPE_STRUCTURES) ;
DATAROM_glo = cell(1,NTYPE) ;
DATA_REFMESH_glo = cell(1,NTYPE) ;
DATAIN.NAME_WS_MODES=cell(1,NTYPE) ;
for itypeLOC = 1:NTYPE
    itype = TYPE_STRUCTURES(itypeLOC) ;
    SCRIPT_MODES = MESH1D.PROP(itype).NameWSmodes;
    
    DATAIN.NAME_WS_MODES{itypeLOC} = [FOLDER,filesep,'MODES',filesep,'MODES_',SCRIPT_MODES,'.mat'];
    cd(EXECUTABLE_FOLDER)
    % LOAD_FROM_MEMORY =1 ;
    % if LOAD_FROM_MEMORY ==1
    disp(['Structural entity = ',num2str(itypeLOC)])
    disp('REtrieving OFFLINE data (DATAROM) ...')
    tic
    load(DATAIN.NAME_WS_MODES{itypeLOC},'DATAROM','DATA_REFMESH') ;
    DATAROM_glo{itype} = DATAROM ;
    DATA_REFMESH_glo{itype} = DATA_REFMESH ;
    toc
    disp('DONE')
end

% --------------------------------------------------------------------
% 1D MESH
% ---------------------------------------------------------------------

% -----------------------------------------------------------
DATAIN.NAME_INPUT_DATA = NAME_INPUT_DATA ;



% STEP 1
% ------
% Assembly stiffness matrix
disp('-------------------------------------')
disp(['Assembly 1D stiffness matrix...'])
tic
[K,ndimINTF ]= AssemblyKbeam(DATAROM_glo,MESH1D,DATAIN) ;
disp(['...Done'])
toc
disp('-------------------------------------')
% -------------------------------------
% Assembly vector of external forces
% -----------------------------------
% INTERFACE FORCES
disp('-------------------------------------')
disp(['Assembly 1D external force vector...'])
tic
P = AssemblyPinterfaces(DATAROM_glo,MESH1D,DATAIN,FORCES,DATA_REFMESH_glo,ndimINTF) ;
% Body and traction forces over elements
[Fdom, fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    	AssemblyFdomains(DATAROM_glo,MESH1D,DATAIN,FORCES,ndimINTF,DATA_REFMESH_glo) ;

F = Fdom + P ;
disp(['Done...'])
disp('-------------------------------------')
toc
% -----------------------------------------------------
% Dirichlet boundary conditions
% ------------------------------
[DOFr,DOFl,aR] = DirichletBNDCondBeam(DATAROM_glo,MESH1D,DISP,DATAIN,ndimINTF) ;



% SOLVING tHE SYSTEM OF reduced EQUATIONS  
disp('-------------------------------------')
disp(['Solving 1D system...'])
tic
a = zeros(size(K,1),1) ;
a(DOFr) = aR ;
a(DOFl) = K(DOFl,DOFl)\(F(DOFl)-K(DOFl,DOFr)*aR) ;
toc
disp(['...Done'])


DATAIN = DefaultField(DATAIN,'SaveInformationROM',[])  ; 
if ~isempty(DATAIN.SaveInformationROM) ; 
    save(DATAIN.SaveInformationROM)  
end

% REACTIONS
Reactions = zeros(size(a)) ; 
Reactions(DOFr) = F(DOFr) - K(DOFr,:)*a ; 

% Plotting interface variables
PlotInterfaceVariablesMOD(DATAIN,a,ndimINTF,MESH1D) ; 

% 3D - RECONSTRUCTION PROCESS
% --------------------------------
% Amplitude self-equilbrated reaction modes
[rDEF]= AmplitudeReactions_jbeam(DATAROM_glo,MESH1D,a,fextBEAMr,ndimINTF) ;

% Axial, shear forces, torsion and bending moments
[GeneralizedForces]= ForcesMoments_diagrams(DATAROM_glo,MESH1D,rDEF,rRB,DATA_REFMESH_glo,ndimINTF,DATAIN) ;

% Amplitude displacement modes
[qDEF,qRB]= AmplitudeDisplacements(DATAROM_glo,MESH1D,rDEF,fextDOMred,DATA_REFMESH_glo,a,DATAIN,ndimINTF) ;

% Error (jumps in displacements). 15th-Aug-2019
% ---------------------------------------------
[DIMENSIONLESS_ERROR]= ErrorInDisplacementsROM...
    (DATAROM_glo,MESH1D,qDEF,qRB,DATA_REFMESH_glo,DATAIN,ndimINTF,DOFr,a) ;





% Reconstruction of displacement and stress fields
% ------------------------------------
disp('-------------------------------------')
disp(['Reconstruction of 3D displacement and stresses...'])
tic
[DISP3D,DISP3D_lateral,STRESS3D,STRESSDATA,DATAIN,DATA_REFMESH_glo,REACTIONS3D]= ...
    Displacement_stress_3D_JOINTslice(DATAROM_glo,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH_glo,...
    DATARUN,rDEF,rRB) ;

toc
disp(['...Done'])



disp('-------------------------------------')
disp(['Printing results in GID ....'])
tic
% Printing post-process file (GID)

% Old impplementation ----just slices. Before July-5th-2018
% GIDprint_BEAM_ROM(MESH1D,DATA_REFMESH,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,STRESS3D,STRESSDATA,...
%     FORCES_2_PRINT) ;
DATAIN.NAME_LOAD_DATA = NAME_LOAD_DATA_inp;
DATAIN.DIMENSIONLESS_ERROR = DIMENSIONLESS_ERROR ; 
GIDprint_BEAM_ROM_JOINTslice(MESH1D,DATA_REFMESH_glo,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,...
    STRESS3D,STRESSDATA,FORCES_2_PRINT,MESH3D,FOLDER,DATARUN,REACTIONS3D,DATAROM_glo)


disp('-------------------------------------')
toc
disp(['...Done'])












