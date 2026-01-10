function ROMnonlinear_straightbeam(NAME_INPUT_DATA,FOLDER,EXECUTABLE_FOLDER,NAME_LOAD_DATA,MESH1D,MESH3D,DATARUN)
if nargin == 0
    load('tmp1.mat')
end
% --------------------------------------
% Reduced-order model for a beam-like structure
% made by repeating 3D slices
% Nonlinear regime. Copy of  ROM_straightbeam
% -------------------------------------------
% JAHO, 9-Oct-2018
% ----------------------------------------------
% --------------------------------------------------------
%  INITIAL  OPERATIONS 
% Data Struture: DATAROM_glo and DATA_REFMESH_glo
% ----------------------------------------------- 
[DATAROM_glo,DATA_REFMESH_glo,DATAIN,FORCES,ndimINTF,MESH1D] ...
    = PrepareData1DbeamROM(MESH1D,FOLDER,NAME_LOAD_DATA,EXECUTABLE_FOLDER,NAME_INPUT_DATA) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EXTERNAL ACTIONS  (ALONG TIME )
% -----------------------------------
% INTERFACE FORCES
disp('-------------------------------------')
disp(['Assembly 1D external force vector...'])
tic
P = AssemblyPinterfacesNON(DATAROM_glo,MESH1D,DATAIN,FORCES,DATA_REFMESH_glo,ndimINTF) ;
% Body and traction forces over elements
[Fdom, fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= AssemblyFdomains(DATAROM_glo,MESH1D,DATAIN,FORCES,ndimINTF) ;
F = Fdom + P ;
disp(['Done...'])




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

% -----------------------------------------

 

toc
% -----------------------------------------------------

% Dirichlet boundary conditions
% ------------------------------
[DOFr,DOFl,aR] = DirichletBNDCondBeam(DATAROM_glo,MESH1D,DISP,DATAIN,ndimINTF) ;

% SOLVING SYSTEM OF EQUATIONS
disp('-------------------------------------')
disp(['Solving 1D system...'])
tic
a = zeros(size(K,1),1) ;
a(DOFr) = aR ;

a(DOFl) = K(DOFl,DOFl)\(F(DOFl)-K(DOFl,DOFr)*aR) ;
toc
disp(['...Done'])

% 3D - RECONSTRUCTION PROCESS
% --------------------------------

% Amplitude self-equilbrated reaction modes
[rDEF]= AmplitudeReactions(DATAROM_glo,MESH1D,a,fextBEAMr,ndimINTF) ;

% Axial, shear forces, torsion and bending moments
[GeneralizedForces]= ForcesMoments_diagrams(DATAROM_glo,MESH1D,rDEF,rRB,DATA_REFMESH_glo,ndimINTF) ;

% Amplitude displacement modes
[qDEF,qRB]= AmplitudeDisplacements(DATAROM_glo,MESH1D,rDEF,fextDOMred,DATA_REFMESH_glo,a,DATAIN,ndimINTF) ;

% Reconstruction of displacement and stress fields
% ------------------------------------
disp('-------------------------------------')
disp(['Reconstruction of 3D displacement and stresses...'])
tic


% DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER = 10 ;
% DATAIN.DOMAINS_POSTPROCESS_SELECT.VARIABLE = 'VONMISES'; % Number of domains to be postprocess




[DISP3D,DISP3D_lateral,STRESS3D,STRESSDATA,DATAIN,DATA_REFMESH_glo]= ...
    Displacement_stress_3D_JOINTslice(DATAROM_glo,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH_glo,DATARUN) ;

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
GIDprint_BEAM_ROM_JOINTslice(MESH1D,DATA_REFMESH_glo,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,...
    STRESS3D,STRESSDATA,FORCES_2_PRINT,MESH3D,FOLDER,DATARUN)


disp('-------------------------------------')
toc
disp(['...Done'])












