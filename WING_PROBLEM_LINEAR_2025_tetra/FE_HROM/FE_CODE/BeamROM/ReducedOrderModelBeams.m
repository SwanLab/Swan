% function ReducedOrderModelBeams(MESH1D,DATAROM,DATA_REFMESH,DATAIN,FORCES,DISP) 
% 
% if nargin == 0
%     load('tmp.mat')
% end
% 
% % STEP 1
% % ------
% % Assembly stiffness matrix
% disp('-------------------------------------')
% disp(['Assembly 1D stiffness matrix...'])
% tic
% K = AssemblyKbeam(DATAROM,MESH1D,DATAIN) ; 
% disp(['...Done'])
% toc
% 
% 
% disp('-------------------------------------')
% 
% % -----------------------------------------
% 
% % -------------------------------------
% % Assembly vector of external forces
% % -----------------------------------
% % INTERFACE FORCES   
% disp('-------------------------------------')
% disp(['Assembly 1D external force vector...'])
% tic
% P = AssemblyPinterfaces(DATAROM,MESH1D,DATAIN,FORCES) ; 
% % Body and traction forces over elements 
% [Fdom, fextBEAMr,rRB,fextDOMred]= AssemblyFdomains(DATAROM,MESH1D,DATAIN,FORCES) ; 
% F = Fdom + P ; 
% disp(['Done...'])
% disp('-------------------------------------')
% 
% toc
% % -----------------------------------------------------
% 
% % Dirichlet boundary conditions 
% % ------------------------------
% [DOFr,DOFl,aR] = DirichletBNDCondBeam(DATAROM,MESH1D,DISP) ;
% 
% % SOLVING SYSTEM OF EQUATIONS   
% disp('-------------------------------------')
% disp(['Solving 1D system...'])
% tic
% a = zeros(size(K,1),1) ; 
% a(DOFr) = aR ; 
% a(DOFl) = K(DOFl,DOFl)\(F(DOFl)-K(DOFl,DOFr)*aR) ; 
% toc
% disp(['...Done'])
% 
% % 3D - RECONSTRUCTION PROCESS 
% % --------------------------------
% 
% % Amplitude self-equilbrated reaction modes
% [rDEF]= AmplitudeReactions(DATAROM,MESH1D,a,fextBEAMr) ;
% 
% % Axial, shear forces, torsion and bending moments 
% [GeneralizedForces]= ForcesMoments_diagrams(DATAROM,MESH1D,rDEF,rRB,DATA_REFMESH) ;
% 
% % Amplitude displacement modes
% [qDEF,qRB]= AmplitudeDisplacements(DATAROM,MESH1D,rDEF,fextDOMred,DATA_REFMESH,a,DATAIN) ; 
% 
% % Reconstruction of displacement and stress fields 
% % ------------------------------------
% disp('-------------------------------------')
% disp(['Reconstruction of 3D displacement and stresses...'])
% tic
% DATAIN = DefaultField(DATAIN,'DOMAINS_POSTPROCESS',1:length(rDEF)) ;
% if isempty(DATAIN.DOMAINS_POSTPROCESS)
%     DATAIN.DOMAINS_POSTPROCESS = 1:length(rDEF) ; 
% end
% toc 
% disp(['...Done'])
% 
%  
% [DISP3D,DISP3D_lateral,STRESS3D,STRESSDATA]= Displacement_stress_3D(DATAROM,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH) ;  
% 
%  disp('-------------------------------------')
% disp(['Printing results in GID ....'])
% tic
% % Printing post-process file (GID)
% GIDprint_BEAM_ROM(MESH1D,DATA_REFMESH,GeneralizedForces,DATAIN,DISP3D,a,DISP3D_lateral,STRESS3D,STRESSDATA) ; 
%  disp('-------------------------------------')
% toc
% disp(['...Done'])
% 
% 
% 
% end
% 
% 
% 
