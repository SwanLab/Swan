function    [DATAIN,NODES_SNAP,DOFsKEEP,NODES_LINES_ROTATIONS,DISP3D]  = ReducedOrderModel_RVE(DATA_REFMESH_glo,DATAROM_glo,MESH2D,DATAIN,FORCES,DISP,MESH3D,FOLDER,DATARUN)

if nargin == 0
    load('tmp3.mat')
end

NODES_SNAP = [] ; 

 [ndimINTF,DATAIN ]= InformationReducedOrderModel(DATAROM_glo,MESH2D,DATAIN) ; 
 DATAIN.ndimINTF  = ndimINTF; 
% STEP 1
% ------
% Assembly stiffness matrix
disp('-------------------------------------')
disp(['Assembly reduced-order stiffness matrix...'])
tic

[K,DOFsKEEP ]= AssemblyKrve(DATAROM_glo,MESH2D,DATAIN,ndimINTF) ;
ndof = size(K,1) ; 

disp(['...Done'])
toc

disp('-------------------------------------')
% -----------------------------------------

% -------------------------------------
% Assembly vector of external forces
% -----------------------------------
% INTERFACE FORCES
disp('-------------------------------------')
disp(['Assembly reduced-order external force vector...'])
tic

P = AssemblyPinterfacesRVE(DATAROM_glo,MESH2D,DATAIN,FORCES,DATA_REFMESH_glo,ndimINTF,DOFsKEEP) ;

% -----------------------------------------
% Body and traction forces over elements
% -----------------------------------------
[Fdom, fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    AssemblyFdomainsRVE(DATAROM_glo,MESH2D,DATAIN,FORCES,ndimINTF,DOFsKEEP,DATA_REFMESH_glo) ;

snapnow;


F = Fdom + P ;
disp(['Done...'])
disp('-------------------------------------')

toc
% -----------------------------------------------------

% Dirichlet boundary conditions
% ------------------------------
DATAIN.ndimSP = size(DATA_REFMESH_glo{1}.COOR,2) ; 
[DOFr,aR,~,~,NODES_LINES_ROTATIONS] = DirichletBNDCondRVE(DATAROM_glo,MESH2D,DISP,DATAIN,ndimINTF,DOFsKEEP,DATA_REFMESH_glo{1}) ;


% SOLVING SYSTEM OF EQUATIONS
disp('-------------------------------------')
tic
DOFl = 1:size(K,1) ; 
disp(['Solving reduced-order equations  (neq = ',num2str(length(DOFl)),')'])

DOFl(DOFr) = [] ; 
a = zeros(size(K,1),1) ;
a(DOFr) = aR ;

a(DOFl) = K(DOFl,DOFl)\(F(DOFl)-K(DOFl,DOFr)*aR) ;
toc
disp(['...Done'])

NODES_SNAP.U = a ; 

% Reactions  
NODES_SNAP.Reactions = zeros(size(a)) ; 
NODES_SNAP.Reactions(DOFr) = K(DOFr,:)*a - F(DOFr) ; 

% 3D - RECONSTRUCTION PROCESS
% --------------------------------

% Amplitude self-equilbrated reaction modes
  
[rDEF]= AmplitudeReactions_RVE(DATAROM_glo,MESH2D,a,fextRVEr,ndimINTF,DOFsKEEP) ;
 

% Axial, shear forces, torsion and bending moments
GeneralizedForces= ForcesMoments_RVE(DATAROM_glo,MESH2D,rDEF,rRB,DATA_REFMESH_glo) ;
% Amplitude displacement modes
[qDEF,qRB]= AmplitudeDisplacements_RVE(DATAROM_glo,MESH2D,rDEF,fextDOMred,...
    DATA_REFMESH_glo,a,DATAIN,ndimINTF,DOFsKEEP) ;
% Reconstruction of displacement and stress fields
% ------------------------------------
disp('-------------------------------------')
disp(['Reconstruction of 3D displacement and stresses...'])
tic
DATAadd = [] ; 
[DISP3D,DISP3D_lateral,STRESS3D,STRESSDATA,DATAIN,DATA_REFMESH_glo]= ...
    Displacement_stress_3D_JOINTrve(DATAROM_glo,MESH2D,qDEF,qRB,DATAIN,DATA_REFMESH_glo,DATAadd,DATARUN) ;
toc
disp(['...Done'])



disp('-------------------------------------')
disp(['Printing results in GID ....'])
tic
% Printing post-process file (GID)

% Old impplementation ----just slices. Before July-5th-2018
% GIDprint_RVE_ROM(MESH2D,DATA_REFMESH,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,STRESS3D,STRESSDATA,...
%     FORCES_2_PRINT) ;

GIDprint_RVE_ROM(MESH2D,DATA_REFMESH_glo,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,...
    STRESS3D,STRESSDATA,FORCES_2_PRINT,MESH3D,FOLDER,DATARUN,DOFsKEEP,ndimINTF,DATAROM_glo)


disp('-------------------------------------')
toc
disp(['...Done'])

