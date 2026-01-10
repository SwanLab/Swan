function[d, stressGLO,    K, M,Ftrac,CN,TRANSF_COORD,DATA,Vrot,stressesREF]  = SolveElastFE_EIFE(COOR,CN,TypeElement,TypeElementB,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,DATA,NameFileMesh,PROPMAT,MaterialType)
% EMPIRICAL INTERSCALE FINITE ELEMENT METHOD
% JAHO, 9-MARCH-2023

%%% This function returns the (nnode*ndim x 1) vector of nodal displacements (d),
%%% as well as the arrays containing  the stresses (stressGLO) and strains (strainGLO) at all gauss
%%% points
% % INPUTS
% --------------
% 1. Finite element mesh
% -------------------
% COOR: Coordinate matrix (nnode x ndim)
% CN: Connectivity matrix (nelem x nnodeE)
% TypeElement: Type of finite element (quadrilateral,...)
% TypeElementB: Type of boundary finite element (linear...)
% -----------
% 2. Material
% -----------
%  celasglo (nstrain x nstrain x nelem)  % Array of elasticity matrices
%  celasgloINV (6 x 6 x nelem)  % Array of compliance matrices (3D)
% -------------------------
% 3. Dirichlet (Essential) Boundary Condition s
% --------------------------------------------
%  DOFr --> Set of Global Degrees of Freedom with prescribed displacements
%  dR   --> Vector of prescribed displacements  (size(DOFr) = size(dR))
% ---------------------------------------
% 4. Neumann (natural) Boundary Conditions
% -------------------------------------------
% DISTRIBUTED LOADS
% -------------------
%  CNb: Cell array in which the cell CNb{idim} contains the connectivity matrix for the boundary elements
%  of the traction boundaries in the idim direction
%  Tnod: Cell array in which the entry Tnod{idim} features the vector with the prescribed traction at
%   the nodes specified in CNb{idim}    (note that size(CNb{idim}) = size(Tnod{idim}))
%  each cell of
%  POINT LOADS
% --------------------------------
%  Fpnt  (nnode*ndime x 1):  Vector containing point forces applied on the
%  nodes of the discretization
% ----------------------------------
% 5. Body force
% ---------------
%  fNOD: Vector containing the nodal values of the heat source function (nnode*ndime x1 )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=[]; strainGLO=[] ; stressGLO=[] ;posgp=[] ;

% ----------------------------------------------------------------------------------------------
% DETERMINATION COORDINATE TRANSFORMATION, JACOBIAN, B-MATRICES,
% N-MATRICES
[Bmat,Nmat,WEIGHTSinteg,TRANSF_COORD,CN,DATA,Vrot  ]= B_N_matricesEIFE(COOR,CN, PROPMAT,MaterialType,DATA,TypeElement) ;
% ---------------------------------------------------------------------------------------------

% COARSE/FINE SCALE MESHES (check conformity all elements)
% -------------------------------------------------------------
DATA = DefaultField(DATA,'CHECK_FINE_MESH',0) ;
DATA.NameFileMesh = NameFileMesh ;
if DATA.CHECK_FINE_MESH == 1
    
    CheckFineScaleMeshGID(COOR,CN, PROPMAT,MaterialType,TypeElement,DATA,TRANSF_COORD) ;
    error('Set DATA.CHECK_FINE_MESH = 0 to continue with the calculations')
    %d= [] ; strainGLO = [] ;  stressGLO = [] ;   React = [] ;  posgp = [] ;  K = [] ;  M = [] ; Ftrac =[] ;
    %return
end


% A) Global stiffness matrix
% ------------------------------
disp('Computing stiffness matrix K ...')
K= ComputeK_EIFEnv(COOR,CN,Bmat,WEIGHTSinteg.INTforces, PROPMAT,MaterialType,TRANSF_COORD) ;

% B) External force vector due to body forces
% ------------------------------
%disp('Computing   external force vector due to body forces (Fb)...')
Fb = 0 ; % ComputeFb(COOR,CN,TypeElement, fNOD);

% C)  External force vector due to   boundary tractions
% ------------------------------
disp('Computing  external force vector due to   boundary tractions ..')
Ftrac = FtracCOMP(COOR,CNb,TypeElementB,Fpnt,Tnod);

% % Thermal problems
% DATA = DefaultField(DATA,'DeltaTemperature',[]) ;
% if ~isempty(DATA.DeltaTemperature)
%     if length(DATA.DeltaTemperature) ~=size(COOR,1)
%         error('Distinct number of temperature entries than nodes')
%     end
%     
%     Fthermal = ComputeFthermal(COOR,CN,TypeElement,DATA.beta_thermalEXP,DATA.DeltaTemperature);
% else
%     Fthermal = 0 ;
% end

%Fb = Fb + Fthermal;

% D) Mass matrix
DATA = DefaultField(DATA,'COMPUTE_MASS_MATRIX',0) ;
if DATA.COMPUTE_MASS_MATRIX == 1
    disp('Computing mass matrix M ...')
    %    M = ComputeM(COOR,CN,TypeElement, densglo) ;
    
    M = ComputeM_EIFEnv(COOR,CN,Nmat,WEIGHTSinteg.BodyForces, PROPMAT,MaterialType,TRANSF_COORD) ;
else
    M = [] ;
end
% D) Solving for the vector of unknown displacements
disp('Solving...')
[d,stressGLO,TRANSF_COORD,stressesREF  ] = ...
    SolveELAS_EIFE(K,Fb,Ftrac,dR,DOFr,COOR,CN,TypeElement,DATA,M,NameFileMesh,PROPMAT,MaterialType,TRANSF_COORD,Bmat,WEIGHTSinteg) ;

