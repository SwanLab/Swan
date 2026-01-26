function [COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType] = INPUT_LAMINATES(INPUTS_LOC,DATA) ;

% Inputs RVE laminate
% ----------------------------
% MACROSCOPIC GENERALIZED STRAINS
% -------------------
% -------------------
%  strainG = zeros(8,1) ;
%     strainG(1) = 0.1;   % ex
%     strainG(2) = 0 ;   % ey
%     strainG(3) = 0;   % ez
%     strainG(4) = 0 ;   %  Curvature k_x
%     strainG(5) = 0 ;   % Curvature k_y
%     strainG(6) = 0; %  Curvature k_xy
%     straingG(7) = 0 ; % shear gamma_yz
%     straingG(8) = 0 ; % shear gamma_yz

if nargin == 0
    load('tmp.mat')
end


strainINP = INPUTS_LOC.strainINP ;

%MACRODEF = Voig2MatSTRAIN(strainINP) ;  % Voig notation

% SOLVER
DATA.TYPESOLVER = 1; % Conjugated gradient
% DATA.niterCONJG = 1000 ; % Number of iterations
% DATA.tolCONJG = 1e-6 ; % Tolerance solver

%%%%%%%%%%%%%%%%%
% ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB


NameFileMesh = INPUTS_LOC.NameFileMesh ;
READ_MATERIAL_COLUMN = 1;


DATA.nameWORKSPACE = ['DATAWS/',NameFileMesh,'_WS.mat'] ; % To retrieve already read data

%dbstop('41')
if exist(DATA.nameWORKSPACE)~=2
    DATA.RECALCULATE_STIFFNESS = 1 ;
end

if DATA.RECALCULATE_STIFFNESS == 1
    disp('Reading GID mesh ...')
    
    [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
        ReadMeshFile(NameFileMesh,'READ_MATERIAL_COLUMN',READ_MATERIAL_COLUMN)  ;
    save(DATA.nameWORKSPACE,'COOR','CN','TypeElement','CONNECTb','TypeElementB','MaterialType')
    
    disp('Done')
else
    disp('Retrieving mesh data... DONE')
    load(DATA.nameWORKSPACE,'COOR','CN','TypeElement','CONNECTb','TypeElementB','MaterialType') ;
end

nnode = size(COOR,1) ;% Number of nodes
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements

%%% RECALCULATE COORDINATES SO THAT THE MIDPLANE --> Z= 0
zmin = min(COOR(:,3)) ; zmin = zmin(1) ;
zmax = max(COOR(:,3)) ; zmax = zmax(1) ;
zmed = 0.5*(zmin+zmax) ;
COOR(:,3) = COOR(:,3)-zmed ;


%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. MATERIAL PROPERTIES: output celasglo   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL.FIBER.Eyoung = 80000 ; %MN/m2
% MATERIAL.FIBER.POISSON = 0.25 ;
% MATERIAL.FIBER.INDEX = 1;
% %%%
% MATERIAL.MATRIX.Eyoung = 20000    ; 200e7 ; %N/m2
% MATERIAL.MATRIX.POISSON = 0.25 ;
% MATERIAL.MATRIX.INDEX = 2;
MATERIAL = INPUTS_LOC.MATERIAL  ;
[celasglo celasgloINV typePROBLEM]  = AssignMatPropLAM(ndim,MATERIAL,nelem,MaterialType) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA = DefaultField(DATA,'BOUNDARY_CONDlam','PERIODIC_uv_TOP_BOTTOM') ;

switch   DATA.BOUNDARY_CONDlam
    case {'PERIODIC'}
        [DOFr DOFm Gb dR dispMACRO NODESpl] = Periodic_BoundaryCOND_lam(COOR,CN,strainINP,CONNECTb,DATA) ;
        
    case {'PERIODIC_uv_TOP_BOTTOM'}
        [DOFr DOFm Gb dR dispMACRO NODESpl] = Periodic_BoundaryCOND_lamTB(COOR,CN,strainINP,CONNECTb,DATA);
     case {'ZERO_minTB'}
        [DOFr DOFm Gb dR dispMACRO NODESpl] = ZERO_minTB_cal(COOR,CN,strainINP,CONNECTb,DATA,TypeElementB);
        
    case {'PERIODIC_minTOPB'}
        
        [DOFr DOFm Gb dR dispMACRO NODESpl] = Periodic_BoundaryCOND_lamMIN(COOR,CN,strainINP,CONNECTb,DATA,TypeElementB);
    case {'PERIODIC_minTOPBzero'}
        [DOFr DOFm Gb dR dispMACRO NODESpl] = Periodic_BoundaryCOND_lamMINzero(COOR,CN,strainINP,CONNECTb,DATA,TypeElementB);
        
    case {'ZERO'}
        [DOFr DOFm Gb dR dispMACRO NODESpl] = ZeroCOND_lam(COOR,CN,strainINP,CONNECTb,DATA) ;
        
    otherwise
        error('Option not implemented')
        
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: Tnod, CNb, Fnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POINT LOADS
% -----------
Fpnt =zeros(ndim*nnode,1) ; % There are no point loads
% DISTRIBUTED LOADS
% ------------------------
CNb ={} ; Tnod={} ;
for idim=1:3
    CNb{idim} = [] ;   % In this case, there are no distributed loads in the idim=1 direction
    Tnod{idim} = [] ;
end
% 5. Body force
%%
fNOD =  zeros(nnode*ndim,1) ;

%%% POST-PROCESS  and other options
DATA.PLOT.REACTIONS = 1 ;
if ~isfield(DATA,'CALCULATE_averageSTRESS')
    DATA.CALCULATE_averageSTRESS = 2 ; % Generalized stresses
end
DATA.VECTcode  = 1 ;   % Vectorize code
DATA.BCBformulation = 1 ;  % Formulation B^T C B
DATA.NOVOIDS = 1 ;
DATA.strainINP  =strainINP ;
DATA.plotFLUCT  =1;
DATA.dispMACRO = dispMACRO;
DATA.NODESpl = NODESpl ;