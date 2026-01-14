function [COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType] = INPUT_LAMINATES_GEN(INPUTS_LOC,DATA) ;

% Inputs   laminate
 %dbstop('7')
if nargin == 0
    load('tmp.mat')
end


DOFm = [] ; Gb = [] ; 

% SOLVER
% DATA.niterCONJG = 1000 ; % Number of iterations
% DATA.tolCONJG = 1e-6 ; % Tolerance solver
 %%%%%%%%%%%%%%%%%
% ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB


NameFileMesh = INPUTS_LOC.NameFileMesh ;
READ_MATERIAL_COLUMN = 1;


DATA.nameWORKSPACE = ['DATAWS/',NameFileMesh,'_WS.mat'] ; % To retrieve already read data

%dbstop('31')
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
[celasglo celasgloINV typePROBLEM densGLO]  = AssignMatPropLAM(ndim,MATERIAL,nelem,MaterialType) ;

DATA.densGLO = densGLO ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rnod = {} ; uPRES={}  ;
% 1) List of nodes at which displacement is prescribed  in DIRECTION  i = 1
% All nodes pertaining to plane z = 0
% Boundary nodes
BoundaryNodes= unique(CONNECTb(:)) ;
xmin = min(COOR(BoundaryNodes,1)) ; xmin = xmin(1) ;
idim=1 ;
rnodBASEloc = find(abs(COOR(BoundaryNodes,1)-xmin)<1e-10) ;
rnodBASE = BoundaryNodes(rnodBASEloc) ;
rnod{idim} =rnodBASE;
% Vector of prescribed displacements
displ1 = 0 ;
uPRES{idim} = displ1*ones(size(rnod{idim})) ;
% 2) List of nodes at which displacement is prescribed  in DIRECTION  i = 2
idim = 2;
rnod{idim} =rnodBASE ;
% Vector of prescribed displacements
displ2 = 0 ;
uPRES{idim} = displ2*ones(size(rnod{idim})) ;
% 2) List of nodes at which displacement is prescribed  in DIRECTION  i = 3
idim = 3;
rnod{idim} =rnodBASE ;
% Vector of prescribed displacements
displ3 = 0 ;
uPRES{idim} = displ3*ones(size(rnod{idim})) ;
%%%% Set of restricted degrees of freedom and vector of prescribed
%%%% displacements (dR)
DOFr = [] ; dR = [] ;
for idim = 1:ndim
    DOFr = [DOFr ; (rnod{idim}-1)*ndim+idim];
    dR = [dR ; uPRES{idim}];
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: Tnod, CNb, Fnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POINT LOADS
% -----------
%dbstop('118')
Fpnt =zeros(ndim*nnode,1) ; % There are no point loads
% DISTRIBUTED LOADS
% ------------------------
CNb =cell(3,1); Tnod=cell(3,1) ;
eval(INPUTS_LOC.NEUMANN_BOUNDARY_CONDITIONS);  
 
%%%%%%%%%%%%%%%%%%%%%%5
 fbody = 0 ; 
fNOD = fbody*ones(nnode*ndim,1) ;


%%% POST-PROCESS  and other options
DATA.PLOT.REACTIONS = 1 ;
     DATA.CALCULATE_averageSTRESS = 0 ; % Generalized stresses
 DATA.VECTcode  = 1 ;   % Vectorize code
DATA.BCBformulation = 1 ;  % Formulation B^T C B
DATA.NOVOIDS = 1 ;
DATA.plotFLUCT  =0;
%DATA.dispMACRO = dispMACRO;
%DATA.NODESpl = NODESpl ;