function [COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType] = INPUTS_GENERAL2D(INPUTS_LOC,DATA) ;

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


[dummy1 NameFileMeshHERE dummy2]= fileparts(NameFileMesh) ; 


DATA = DefaultField(DATA,'INPUTDATAfile',NameFileMeshHERE)  ; 
DATA.nameWORKSPACE = ['DATAWS/',DATA.INPUTDATAfile,'_WS.mat'] ; % To retrieve already read data

%dbstop('31')
if exist(DATA.nameWORKSPACE)~=2
    DATA.RECALCULATE_STIFFNESS = 1 ;
end

DATA= DefaultField(DATA,'FactorDivideCoordinates',1) ; 




if DATA.RECALCULATE_STIFFNESS == 1
    disp('Reading GID mesh ...')    
    
    DATA = DefaultField(DATA,'MakeMeshByRepetition',[]) ; 
    DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'nDOMx',[]) ; 

    if isempty(DATA.MakeMeshByRepetition.nDOMx)
    [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
        ReadMeshFile(NameFileMesh,'READ_MATERIAL_COLUMN',READ_MATERIAL_COLUMN)  ;
    
    else
        nDOM = DATA.MakeMeshByRepetition.nDOMx ; 
         [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
             MeshGenerationRepeat(NameFileMesh,nDOM,DATA);
    end
    
    
    
    
    COOR = COOR/DATA.FactorDivideCoordinates ; 
    COOR = COOR(:,1:2) ; 
    save(DATA.nameWORKSPACE,'COOR','CN','TypeElement','CONNECTb','TypeElementB','MaterialType','NameFileMesh')
        disp('Done')
else
    disp('Retrieving mesh data... DONE')
    load(DATA.nameWORKSPACE,'COOR','CN','TypeElement','CONNECTb','TypeElementB','MaterialType') ;
end







nnode = size(COOR,1) ;% Number of nodes
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements
 


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
%dbstop('81')
eval(INPUTS_LOC.DISPLACEMENT_COND_FILE);



%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: Tnod, CNb, Fnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS and POINT LOADS
% ------------------------
%dbstop('122')
eval(INPUTS_LOC.NEUMANN_BOUNDARY_CONDITIONS);

%%%%%%%%%%%%%%%%%%%%%%5

% Body forces 
% -------------
DATA = DefaultField(DATA,'INCLUDE_SELFWEIGHT',0) ; 
if DATA.INCLUDE_SELFWEIGHT == 1
    fNOD = zeros(nnode*ndim,1) ;
    fNOD(2:ndim:end) = DATA.selfweight  ; 

else
    fbody = 0 ;
fNOD = fbody*ones(nnode*ndim,1) ;
end
 
  
 

%fNOD(933*3) = 0.1 ; % Mpa


%%% POST-PROCESS  and other options
DATA.PLOT.REACTIONS = 1 ;
DATA.CALCULATE_averageSTRESS = 0 ; % Generalized stresses
DATA.VECTcode  = 1 ;   % Vectorize code
DATA.BCBformulation = 1 ;  % Formulation B^T C B
DATA.NOVOIDS = 1 ;
DATA.plotFLUCT  =0;
%DATA.dispMACRO = dispMACRO;
%DATA.NODESpl = NODESpl ;