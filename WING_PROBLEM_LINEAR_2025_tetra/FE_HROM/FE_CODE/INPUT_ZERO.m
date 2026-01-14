function [COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType] = INPUT_ZERO(INPUTS_LOC,DATA) ;

% Inputs UNIT CELL
% ----------------------------
% MACROSCOPIC STRAINS
% -------------------
% ------------------- voig's notation
% strainINP = zeros(6,1) ;
% strainINP(1) = 0;   % ex
% strainINP(2) = 0 ;   % ey
% strainINP(3) = 1e-3 ;   % ez
% strainINP(4) = 0 ;   % gamma_yz
% strainINP(5) = 0 ;   % gamma_xz
% strainINP(6) = 0; % gamma_xy

if nargin == 0
    load('tmp.mat')
end


strainINP = INPUTS_LOC.strainINP ;

MACRODEF = Voig2MatSTRAIN(strainINP) ;  % Voig notation

% SOLVER
%DATA.TYPESOLVER = 1; % Conjugated gradient
DATA = DefaultField(DATA,'TYPESOLVER',1) ;
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
[celasglo celasgloINV typePROBLEM]  = AssignMatProp_ElastIso(ndim,MATERIAL,nelem,MaterialType) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DOFr DOFm Gb dR dispMACRO] = ZERO_KBCS_cube(COOR,CN,MACRODEF,CONNECTb,DATA) ;

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
DATA = DefaultField(DATA,'PLOT',[]) ; 
DATA.PLOT = DefaultField(DATA.PLOT,'REACTIONS',0) ; 

DATA.CALCULATE_averageSTRESS = 1 ;
DATA.VECTcode  = 1 ;   % Vectorize code
DATA.BCBformulation = 1 ;  % Formulation B^T C B
DATA.NOVOIDS = 1 ;
DATA.strainINP  =strainINP ;
DATA.plotFLUCT  =1;
DATA.dispMACRO = dispMACRO;