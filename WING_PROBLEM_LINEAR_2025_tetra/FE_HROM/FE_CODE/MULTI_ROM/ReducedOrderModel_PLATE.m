function ReducedOrderModel_PLATE(DATA_REFMESH_glo,DATAROM_glo,MESH2D,DATAIN,FORCES,DISP,MESH3D,FOLDER,DATARUN)




% STEP 1
% ------
% Assembly stiffness matrix
disp('-------------------------------------')
disp(['Assembly 2D stiffness matrix...'])
tic
[K,ndimINTF ]= AssemblyKrve_PLATE(DATAROM_glo,MESH2D,DATAIN,DATA_REFMESH_glo) ;
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
P = AssemblyPinterfacesRVE_PLATE(DATAROM_glo,MESH2D,DATAIN,FORCES,DATA_REFMESH_glo,ndimINTF) ;
% -----------------------------------------
% Body and traction forces over elements
% -----------------------------------------
[Fdom, fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    AssemblyFdomainsRVE_PLATE(DATAROM_glo,MESH2D,DATAIN,FORCES,ndimINTF) ;

F = Fdom + P ;
disp(['Done...'])
disp('-------------------------------------')
toc
% -----------------------------------------------------
% Dirichlet boundary conditions
% ------------------------------
[DOFr,DOFl,aR] = DirichletBNDCondRVE_PLATE(DATAROM_glo,MESH2D,DISP,DATAIN,ndimINTF,DATA_REFMESH_glo) ;
% SOLVING SYSTEM OF EQUATIONS
disp('-------------------------------------')
disp(['Solving reduced-order equations...'])
tic
a = zeros(size(K,1),1) ;
a(DOFr) = aR ;
a(DOFl) = K(DOFl,DOFl)\(F(DOFl)-K(DOFl,DOFr)*aR) ;
toc
disp(['...Done'])
% 3D - RECONSTRUCTION PROCESS
% --------------------------------
% Amplitude self-equilbrated reaction modes
[rDEF]= AmplitudeReactions_PLATE(DATAROM_glo,MESH2D,a,fextRVEr,ndimINTF,DATA_REFMESH_glo) ;
% Axial, shear forces, torsion and bending moments
GeneralizedForces= ForcesMoments_PLATE(DATAROM_glo,MESH2D,rDEF,rRB,DATA_REFMESH_glo,ndimINTF) ;
% Amplitude displacement modes
[qDEF,qRB]= AmplitudeDisplacements_PLATE(DATAROM_glo,MESH2D,rDEF,fextDOMred,DATA_REFMESH_glo,a,DATAIN,ndimINTF) ;
% Reconstruction of displacement and stress fields
% ------------------------------------
disp('-------------------------------------')
disp(['Reconstruction of 3D displacement and stresses...'])
tic

[DISP3D,DISP3D_lateral,STRESS3D,STRESSDATA,DATAIN,DATA_REFMESH_glo]= ...
    Displacement_stress_3D_JOINTrve(DATAROM_glo,MESH2D,qDEF,qRB,DATAIN,DATA_REFMESH_glo) ;
toc
disp(['...Done'])



disp('-------------------------------------')
disp(['Printing results in GID ....'])
tic
% Printing post-process file (GID)

% Old impplementation ----just slices. Before July-5th-2018
% GIDprint_RVE_ROM(MESH2D,DATA_REFMESH,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,STRESS3D,STRESSDATA,...
%     FORCES_2_PRINT) ;

GIDprint_PLATE_ROM(MESH2D,DATA_REFMESH_glo,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,...
    STRESS3D,STRESSDATA,FORCES_2_PRINT,MESH3D,FOLDER,DATARUN)


disp('-------------------------------------')
toc
disp(['...Done'])

