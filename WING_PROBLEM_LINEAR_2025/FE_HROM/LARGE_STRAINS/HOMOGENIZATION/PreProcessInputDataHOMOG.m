function   [MESH,MATPRO,OPERFE,DISP_CONDITIONS,DATA,OTHER_data] ...
    =  PreProcessInputDataHOMOG(DATA,PROPMAT,MACRODEF,OTHER_INPUTS)
% See  08_ComputationalHomogenization.pdf
% Inputs for homogenization  problems 
% JAHO, 17-JAN-2021
 

 %
OTHER_data =[] ; 
if nargin == 0
    load('tmp1.mat')
end

 

% ----------------------------
%%%%%%%%%%%%%%%%%
% ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB, Connectivity faces,
% normals ....
MESH = GeometryMesh(DATA) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nnode,ndim ]= size(MESH.COOR) ;% Number of nodes
[nelem,nnodeE ]= size(MESH.CN) ; % Number of elements


DATA.MESH.nnode =nnode;
DATA.MESH.ndim =ndim ;
DATA.MESH.nelem =nelem ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. MATERIAL PROPERTIES: output celasglo, density   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MATPRO] = SmallStrainElasticityPROP(MESH,DATA.typePROBLEM,PROPMAT) ;  

% 3. GEOMETRIC MATRICES (AND RENUMBERING)
% Geometric matrices 
%-------------------------------------
DATA = DefaultField(DATA,'RenumberElementsForEficiency',1) ;  
if  DATA.RenumberElementsForEficiency == 1
    disp('Renumering elements')
    [~,IndicesRenumberingElements]  = sort(MESH.CN(:,1)) ;
    MESH.CN = MESH.CN(IndicesRenumberingElements,:) ;
    MATPRO.celasglo = MATPRO.celasglo(:,:,IndicesRenumberingElements) ;
    MATPRO.dens = MATPRO.dens(IndicesRenumberingElements) ;
        MESH.MaterialType = MESH.MaterialType(IndicesRenumberingElements) ;
else
    IndicesRenumberingElements = [] ;
end
MESH.IndicesRenumberingElements =IndicesRenumberingElements ; 
nstrain = size(MATPRO.celasglo,2) ; 
DATA.MESH.nstrain  =nstrain ; 
MESH.DATA = DATA; 
% 4. FE/HROM OPERATORS,MATRICES 
% ---------------------------
[Bst_F,wSTs,Nst,wSTs_RHS,NstT_W_N_boundaries,ngaus_RHS,GEOproperties,ngaus_STRESS,IDENTITY_F,posgp] ...
    = GeometricMatricesFun(MESH,nstrain) ; 

DATA = DefaultField(DATA,'STORE_MATRIX_SHAPE_FUNCTIONS_Nst',1) ;
if DATA.STORE_MATRIX_SHAPE_FUNCTIONS_Nst == 1
    OTHER_data.Nst = Nst ;
end

OTHER_data.GEOproperties = GEOproperties; 
DATA.MESH.posgp  =posgp ; 

DATA = DefaultField(DATA,'TYPE_CONSTITUTIVE_MODEL_ALL','SMALL_STRAINS_ELASTIC') ;  
DATA.MESH.ngaus  =ngaus_STRESS ; 
switch DATA.TYPE_CONSTITUTIVE_MODEL_ALL
    case 'SMALL_STRAINS_ELASTIC'
        % We can precompute the diagonal elasticity matrix at all the Gauss points
        %-------------------------------------------------
        % Matrix with all elasticity tensors (at all gauss points)
        disp('Global elasticity matrix    ...')
         Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...        
      %   Cglo = ConvertBlockDiag(Cglo) ;
      %   MATPRO.celasglo = [] ; 
         MATPRO.celasglo = Cglo ; 
        
end


%OPERFE.Bst = Bst_F ;   % Matrix such that F = Bst_F*d + ident
OPERFE.IDENTITY_F = IDENTITY_F ;
%OPERFE.wSTs = wSTs ;     % Integration weights at the Gauss points used for evaluating stresses and internal
                              % forces 

  OPERFE.wSTs = wSTs ;                
  
% Homogeneization 
% ***************
 
                              
DATA.MESH.ngausT =  length(wSTs) ; % Total number of Gauss points
DATA.MESH.ndofSTRESS =   DATA.MESH.ngausT*nstrain ; % Total number of Gauss points
DATA.MESH.ngaus_STRESS = ngaus_STRESS; 
DATA.MESH.ndof = nnode*ndim; 
DATA.MESH.ndim =  ndim; 

% 5.  MASS MATRIX 
% -----------------------
densGLO = repmat(MATPRO.dens',ngaus_RHS,1) ;
densGLO = densGLO(:) ;
densGLO_W = CompWeightDiag(densGLO.*wSTs_RHS,ndim) ;  
NstW = densGLO_W*Nst ;
%OPERFE.M = NstW'*Nst ;

% COMPUTING NUMBER OF CLUSTERS (FOR STORING PURPOSES)
% ******************************
DATA = DefaultField(DATA,'LIMIT_mbytes_matrices',50) ; 
DATA = ClusterComputeSize(DATA.LIMIT_mbytes_matrices,DATA) ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER_INPUTS = DefaultField(OTHER_INPUTS,'RB_MOTION',[]) ; 
% [DOFr,dR] = DirichletCONDtime(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ; 
% DISP_CONDITIONS.DOFr = DOFr; 
% DOFl = 1:DATA.MESH.ndof ; 
% DOFl(DOFr) = []  ; 
% DISP_CONDITIONS.DOFl = DOFl ;  
% DISP_CONDITIONS.dR = dR; 
% dR = DefaultField(dR,'RIGID_BODY_MOTION',[]) ; 
% DISP_CONDITIONS.RIGID_BODY_MOTION = dR.RIGID_BODY_MOTION ; 
% dR.RIGID_BODY_MOTION = [] ; 

[DISP_CONDITIONS,COORrel,DATA] = Periodic_BoundaryCOND_LARGE(MESH.COOR,MESH.CN,MESH.CNb,DATA) ;
MESH.COOR = COORrel' ; 
MACROVAR = DispMACROtime(MACRODEF,MESH.COOR',DATA) ; 
DISP_CONDITIONS.MACROVAR = MACROVAR ; 

OPERFE.BstA = Bst_F*DISP_CONDITIONS.A ; 

% HOMOGENIZATION OPERATOR 
VOL = DATA.VOL_RVE  ; 
VOL_SOLID = sum(wSTs) ; 
POR = (VOL-VOL_SOLID)/VOL*100 ;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Porososity = ',num2str(POR),' %']) ; 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
OPERFE.WEIGHTShomog = wSTs./VOL ; 
 

 
 
 
% ---------------- 
