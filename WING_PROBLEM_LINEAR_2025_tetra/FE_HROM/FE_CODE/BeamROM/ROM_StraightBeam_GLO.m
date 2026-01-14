function ROM_StraightBeam_GLO(SCRIPT_MODES,EXECUTABLE_FOLDER,DATAIN,FORCES,DISP,NSLICES,NAME_INPUT_DATA)

if nargin == 0
    load('tmp.mat')
end

% Computing stiffness matrix and other reduced-order variables.
% -------------------------------------------------------------
DATAIN.NAME_WS_MODES = [cd,filesep,'MODES',filesep,'MODES_',SCRIPT_MODES,'.mat'];

cd(EXECUTABLE_FOLDER)
LOAD_FROM_MEMORY =1 ;
if LOAD_FROM_MEMORY ==1
    disp('REtrieving OFFLINE data (DATAROM) ...')
    tic
    load(DATAIN.NAME_WS_MODES,'DATAROM','DATA_REFMESH')
    toc
    disp('DONE')
else
    % LOAD BASIS MATRICES
    % -------------------
    nBASES_BEAM = [] ;
    load(DATAIN.NAME_WS_MODES,'BASES','DATA_REFMESH','nBASES_BEAM') ;
    DATAROM =  BeamStiffnessMatrix(BASES,DATA_REFMESH,DATAIN,nBASES_BEAM) ;
end
%
% --------------------------------------------------------------------
% 1D MESH
% ---------------------------------------------------------------------
% Matrix of coordinates for interfaces (1D problem)
L = DATA_REFMESH.LENGTH ;  % Length of the slice in the x-direction
x = 0:L:NSLICES*L ;
MESH1D.COOR = zeros(length(x),3) ;
MESH1D.COOR(:,1) = x ;
MESH1D.LEFT_END_NODE = 1 ; MESH1D.RIGHT_END_NODE = length(x) ;
% Matrix of connectivities
MESH1D.CN = [(1:(length(x)-1))', (2:length(x))'] ;
% MaterialType 1D
MESH1D.MaterialType = ones(size(MESH1D.CN,1),1) ;
% -----------------------------------------------------------
DATAIN.NAME_INPUT_DATA = NAME_INPUT_DATA ;

%ReducedOrderModelBeams(MESH1D,DATAROM,DATA_REFMESH,DATAIN,FORCES,DISP) ;


% STEP 1
% ------
% Assembly stiffness matrix
disp('-------------------------------------')
disp(['Assembly 1D stiffness matrix...'])
tic
K = AssemblyKbeam(DATAROM,MESH1D,DATAIN) ;
disp(['...Done'])
toc


disp('-------------------------------------')

% -----------------------------------------

% -------------------------------------
% Assembly vector of external forces
% -----------------------------------
% INTERFACE FORCES
disp('-------------------------------------')
disp(['Assembly 1D external force vector...'])
tic
P = AssemblyPinterfaces(DATAROM,MESH1D,DATAIN,FORCES,DATA_REFMESH) ;
% Body and traction forces over elements
[Fdom, fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= AssemblyFdomains(DATAROM,MESH1D,DATAIN,FORCES) ;
F = Fdom + P ;
disp(['Done...'])
disp('-------------------------------------')

toc
% -----------------------------------------------------

% Dirichlet boundary conditions
% ------------------------------
[DOFr,DOFl,aR] = DirichletBNDCondBeam(DATAROM,MESH1D,DISP,DATAIN) ;

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
[rDEF]= AmplitudeReactions(DATAROM,MESH1D,a,fextBEAMr) ;

% Axial, shear forces, torsion and bending moments
[GeneralizedForces]= ForcesMoments_diagrams(DATAROM,MESH1D,rDEF,rRB,DATA_REFMESH) ;

% Amplitude displacement modes
[qDEF,qRB]= AmplitudeDisplacements(DATAROM,MESH1D,rDEF,fextDOMred,DATA_REFMESH,a,DATAIN) ;

% Reconstruction of displacement and stress fields
% ------------------------------------
disp('-------------------------------------')
disp(['Reconstruction of 3D displacement and stresses...'])
tic


% DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER = 10 ;
% DATAIN.DOMAINS_POSTPROCESS_SELECT.VARIABLE = 'VONMISES'; % Number of domains to be postprocess




[DISP3D,DISP3D_lateral,STRESS3D,STRESSDATA,DATAIN,DATA_REFMESH]= ...
    Displacement_stress_3D(DATAROM,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH) ;

toc
disp(['...Done'])



disp('-------------------------------------')
disp(['Printing results in GID ....'])
tic
% Printing post-process file (GID)
GIDprint_BEAM_ROM(MESH1D,DATA_REFMESH,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,STRESS3D,STRESSDATA,...
    FORCES_2_PRINT) ;
disp('-------------------------------------')
toc
disp(['...Done'])



end




