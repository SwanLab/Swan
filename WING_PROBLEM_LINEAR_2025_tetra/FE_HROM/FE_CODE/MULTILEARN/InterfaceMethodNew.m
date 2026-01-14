function [maxD] =InterfaceMethodNew(NAMEWS,DATAONLINE,DATAINM)
% Copy of OnlineDomainSVD_INTERFACEm
% Structural analysis using domain-wise reduced-order modeling
% See MultiLearn.pdf, Implementation.pdf
%---------------------------------------------------------------
% INPUTS
%-----------------------------------------------------------
%dbstop('9')
if nargin == 0
    load('tmp1.mat')
end
fff = fieldnames(DATAONLINE) ;
for ifield = 1:length(fff) ;
    STRG = [fff{ifield},'=','DATAONLINE.',fff{ifield},';'] ;
    eval(STRG) ;
end
if exist('INPUT_PERIODIC')==0 ;     addpath('FE_CODE') ; end
if exist(PROJECT_LOADS{1}) ==0 ;     addpath('DATA_input'); end
% ----------
% STEP 1
% ---------
% Retrieving OFFLINE information
load(NAMEWS,'BasisUdef','BasisRdef','BasisUrb','CNref',...
    'COORref','NODESfacesLOC','TypeElement','posgp','Wdom','BdomRED','SingVal_disp','DATAOUT',...
    'Ui_Si','SingVal_reac','Ui_Si_reac') ;
DATAINM = DefaultField(DATAINM,'CUBATURE',[]) ;
DATAINM.CUBATURE = DefaultField(DATAINM.CUBATURE,'ACTIVE',0) ;
 
   
% LOOP OVER PROJECTS DEFINING EXTERNAL FORCES AND STIFFNESS MATRIX
% -----------------------------------------------------------------
%
[KdomRED,FORCE_PROJECTS,DATA,Cglo,MaterialType,nMAT,K] = ExtractPropertiesProjects...
    (PROJECT_LOADS,COORref,BdomRED,Wdom,DATAINM) ;
% -------------------------------------------------------------------------
% STEP 3  - Determination interface nodes and Dirichlet Boundary conditions
% -------------------------------------------------------------------------
ndim  = size(COORref,2) ;
% Determination of nodes pertaining to faces, edges and corners
[NODESbound,DATALIM] =  PointPlanesRBODY_GEN(COORref,CNref,DATA) ;
NODESfaces = NODESbound.PLANE ;

load(NAMEWS,'BasisRrb');
%end
% Dirichlet conditions
uBAR = DomainRedDirichledBC_2D3Dtiling(alphaBC,uBAR,NODESfaces,BasisUrb,ndim,COORref,DATAONLINE,DATALIM) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of RVE for each domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDOM = size(alphaBC,1) ;
[COLUMNS_RVEloc,nTYPErve,COLUMNS_RVE] = SeparatedModesOperations(DATAINM,nDOM) ;
% For later purposes, we define   global basis matrices
% ----------------------------------------------------
[BasisUrbGLO,BasisUdefGLO,BasisRdefGLO,BasisRrbGLO] = GlobalBasisMatrices_INTFM(BasisRrb,BasisUrb,nDOM,BasisUdef,BasisRdef,COLUMNS_RVE);

% STEP 4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Assembly system of equations   A x = b
% % ----------------------------------------
% -----------------------------------
% A) Matrix Pcomp (compatibility matrix ), and vectors bCOMPr% qRB = PcompRR\bCOMPr ;
%  and bCOMPd
% % % ---------------------------------
% %[PcompRR PcompRD PcompDD bCOMPr bCOMPd] =  AssemblyPcompNEW(BasisUrb,BasisUdef,f1,f2,nDOM,uBAR,ALPHA_ENDS) ;
% [Pcomp, bCOMP, mCOMP,INDrig,INDdef DATAwmethod] = ...
%     AssemblyPcompTILING2D3D(BasisUrb,BasisUdef,NODESfaces,uBAR,alphaBC,COLUMNS_RVE,COLUMNS_RVEloc,ndim,...
%     THETAfaces,GAMMAfaces,DATAINM) ;

% B) Rigid body reactions
% % RB reactions
[fextDOMglo,rRBglo,reactDOMrbGLO] = RigidBodyReactions(BasisUrb,nDOM,BasisRrb,FORCE_PROJECTS,INTENSITY_LOADS) ;

%%%%%%%%%%%%%%%%%
% C) Transference matrix
% ----------------------
% % disp('Assembly of Treac')
% [Treac, cREAC, s, INDrigR,INDdefR,Tall]=  AssemblyTreacTILING2D3D(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
%     THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D) Reaction-displacement basis matrix
% -----------------------------
FACES_CONTACT = sum(THETAfaces,1) ;
IND_CONTACT= find(FACES_CONTACT) ;
if ~isempty(IND_CONTACT) && sum(IND_CONTACT-[1 3]) ~=0
    error('Option not compatible. ONly 1D-tiling problems are allowed (along x-axis)')
end
[KdomREDglo,Hdr,Hrd,Hrr,Hcond] = RedStiffnessMatrix_Hmatrix_INTF(BasisUdef...
    ,nDOM,COLUMNS_RVE,KdomRED,BasisRdef,DATAINM,BasisUrb,BasisRrb,alphaBC,NODESbound,ndim) ;

%
% F) External forces
fextDOMglo_mat  = cell2mat(fextDOMglo) ;
fextDOMglo_mat = (fextDOMglo_mat(:)) ;
fextRED_d = BasisUdefGLO'*fextDOMglo_mat  ;
fextRED_r = BasisUrbGLO'*fextDOMglo_mat  ;

% G) Matrix reaction modes -interface displacement modes (Tint,BasisINT).
% Boundary terms
[Tcond,Tr ,cBC_d, cBC_r,BasisINT   ] =   Tint_InterfaceMethod(BasisUdef,nDOM,BasisRdef,...
    DATAINM,BasisUrb,BasisRrb,alphaBC,NODESbound,ndim,COORref,uBAR,Hdr,Hrd,Hrr,...
    CNref,SingVal_disp,DATAINM,DATAOUT,Ui_Si,SingVal_reac,Ui_Si_reac) ;



% H) Solving the reduced equations 
[qDEF,qRB,rDEF,rRB,pINT] = EquationsAssemblyInterfaceM(Tcond,Tr,cBC_d,cBC_r,KdomREDglo,...
    Hdr,Hrd,Hrr,Hcond,fextRED_d,fextRED_r) ;
% I) COmputing residual

[RESIDUAL,DATAINM ]= ComputeResidual(qDEF,rDEF,rRB,fextDOMglo,BasisUdef,BasisRdef,BasisRrb,K,DATAINM);

% Displacement field
d = BasisUrbGLO*qRB + BasisUdefGLO*qDEF ;
disp('MAx  y disp =')
disp(['MAxDISPY =',num2str(max(abs(d))'), ' (m)'])
REACT =  BasisRdefGLO*rDEF + BasisRrbGLO*rRB  ;

%
% end
% Stress reconstruction
if exist('setPoints')==0
    setPoints = [];
end
DATAINM = DefaultField(DATAINM,'PRINT_AVERAGE_STRESSES_ON_ELEMENTS',1) ;
stressVONMISES= [] ; MAXstressVONMISES = [] ;
% [stressGLO,setElementsRED] = ...
%     StressReconstructionMULTI(qDEF,posgp,setPoints,BdomRED,Cglo,nDOM,DATAINM,Wdom,ndim,COLUMNS_RVE,COLUMNS_RVEloc) ;

%%%  PRINTING OPTIONS
%  DEFAULT VALUES
DATAINM = PrintingOptions(DATAINM) ;

% [stressGLO,setElementsRED,stressVONMISES,MAXstressVONMISES] = ...
%     StressReconstructionMULTI_elem(qDEF,posgp,setPoints,BdomRED,Cglo,nDOM,DATAINM,Wdom,ndim,COLUMNS_RVE,COLUMNS_RVEloc) ;
[stressGLO,setElementsRED,stressVONMISES,MAXstressVONMISES,DOMAINS_TO_INCLUDE] = ...
    StressReconstructionMULTI_CELL(qDEF,posgp,setPoints,BdomRED,Cglo,nDOM,DATAINM,Wdom,ndim,COLUMNS_RVE,COLUMNS_RVEloc) ;
DATA.setElementsRED = setElementsRED ;



%%%% GENERATING FINITE ELEMENT MESH by repeating unit cells

[COORrecons,CNrecons,Materials,dFACES,DOFs_to_INCLUDE] = ...
    MeshByRepetitionTILING2D3D(COORref,CNref,NODESfaces,DATAONLINE,MaterialType,DATAINM,nMAT,...
    DOMAINS_TO_INCLUDE) ;


% -----------------------------------------------------------------
%%% PRINTING GID RESULTS
% -----------------------
%if iproj == length(nDEFglo)
nDEF = size(BasisUdef{1},2) ;
NAME_BASE_GIDfiles = [NAME_PROJECT_TEST,'_APPROX','_',num2str(nDEF)] ;
NameFileMeshss = [] ;
ndim =  size(COORrecons,2) ;
DATAINM = DefaultField(DATAINM,'PlotGidTypeDisplacements','total') ;
switch DATAINM.PlotGidTypeDisplacements
    case {'total',''}
    case 'RIGID'
        d = BasisUrbGLO*qRB;
    case 'STRAIN'
        d = BasisUdefGLO*qDEF ;
    otherwise
        error('Option not implemented')
end
if ~isempty(DOFs_to_INCLUDE)
    d = d(DOFs_to_INCLUDE);
    REACT = REACT(DOFs_to_INCLUDE) ;
    if ~isempty(RESIDUAL.ALL)
        RESIDUAL.ALL = RESIDUAL.ALL(DOFs_to_INCLUDE) ; 
    end
end


strainGLO = [] ;
DATAINM = DefaultField(DATAINM,'OPEN_GID',0) ; % = 1;
DATA.OPEN_GID =DATAINM.OPEN_GID ; %  = DefaultField(DATAIN,'OPEN_GID',0) ; % = 1;
DATA.stressVONMISES = stressVONMISES ;
if DATAINM.GID_PRINT_ALL_DOMAINS==1
    DATA.MakeMeshByRepetition.nDOM = DATAONLINE.NdomX ;
    GidPostProcess(COORrecons,CNrecons,TypeElement,d,strainGLO, ...
        stressGLO,  REACT,NAME_BASE_GIDfiles,posgp,NameFileMeshss,Materials,DATA,...
        RESIDUAL);
end
%end
if ~isempty(setElementsRED)
    clipboard('copy',num2str(setElementsRED'));
end

%
dNORM = reshape(d,ndim,[]) ;
dNORM = sqrt(sum(dNORM.*dNORM,1)) ;
maxD = max(dNORM) ;
% PRINT 1D REPRESENTATION
DATAONLINE = DefaultField(DATAONLINE,'NdomZ',1) ;
DATAONLINE = DefaultField(DATAONLINE,'NdomY',1) ;
if DATAINM.PRINT_GID_1D_representation == 1  & DATAONLINE.NdomY ==1 & DATAONLINE.NdomZ ==1
    PrintGID1Drepresentation(DATAINM,dFACES,DATAONLINE,pINT,alphaBC,BasisINT,COORref,NODESfaces,...
        BasisRdef,BasisRrb,rDEF,rRB,MAXstressVONMISES) ;
end
if DATAINM.GID_1D_3D_print.ACTIVE == 1  & DATAONLINE.NdomY ==1 & DATAONLINE.NdomZ ==1
    PrintGID1D_3Drepresentation(DATAINM,dFACES,DATAONLINE,pINT,alphaBC,BasisINT,COORref,NODESfaces,...
        BasisRdef,BasisRrb,rDEF,rRB,MAXstressVONMISES,COORrecons,CNrecons,TypeElement,d,...
        strainGLO,stressGLO,  REACT,NAME_BASE_GIDfiles,posgp,NameFileMeshss,Materials,DATA,RESIDUAL) ;
end
