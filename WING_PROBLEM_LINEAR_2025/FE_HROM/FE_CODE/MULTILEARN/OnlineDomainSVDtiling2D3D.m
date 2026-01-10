function [Q_d,Q_r,maxD] =OnlineDomainSVDtiling2D3D(NAMEWS,DATAONLINE,DATAINM)
% Structural analysis using domain-wise reduced-order modeling
% Prototype for one single domain
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
    'COORref','NODESfacesLOC','TypeElement','posgp','Wdom','BdomRED','DOFS_reference') ;
DATAINM = DefaultField(DATAINM,'CUBATURE',[]) ;
DATAINM.CUBATURE = DefaultField(DATAINM.CUBATURE,'ACTIVE',0) ;
DATAINM = DefaultField(DATAINM,'MinimizationBoundaryWork',0) ;
if DATAINM.MinimizationBoundaryWork == 1
error('DATA.MinimizationBoundaryWork=0 is conducive to ill-posed problems, see Implementation.pdf')
end
DATAINM = DefaultField(DATAINM,'ConstrainedWithReactionsAtInterfaces',0) ;

if DATAINM.CUBATURE.ACTIVE == 1
    load(NAMEWS,'BasisS','BdomRED','setPoints','WdomRED')
else
    setPoints = [] ; WdomRED = [] ; 
end
% ----------
% STEP 2
% ---------
% LOOP OVER PROJECTS DEFINING EXTERNAL FORCES AND STIFFNESS MATRIX
% -----------------------------------------------------------------
%
[KdomRED,FORCE_PROJECTS,DATA,Cglo,MaterialType,nMAT,K] = ExtractPropertiesProjects...
(PROJECT_LOADS,COORref,BdomRED,Wdom,DATAINM,setPoints,WdomRED) ;
% -------------------------------------------------------------------------
% STEP 3  - Determination interface nodes and Dirichlet Boundary conditions
% -------------------------------------------------------------------------
ndim  = size(COORref,2) ;
% Determination of nodes pertaining to faces, edges and corners
[NODESbound,DATALIM] =  PointPlanesRBODY_GEN(COORref,CNref,DATA) ;
NODESfaces = NODESbound.PLANE ;
% DATAIN.RecalculateRigidBodyMatrix = 0;
% if DATAIN.RecalculateRigidBodyMatrix  == 1
%     f = unique(cell2mat(NODESfaces)) ;
%     f = small2large(f,ndim);
%     % Rigid body-reactions mode matrix
%     BasisRrb = zeros(size(BasisUrb)) ;
%     BasisRrb(f,:) = BasisUrb(f,:) ;
% else
% Taking all DOFs by default as restricted surfaces may give rise to
% problems. In this case, it is better to retrieve it from memory
load(NAMEWS,'BasisRrb');
%end
% Dirichlet conditions
uBAR = DomainRedDirichledBC_2D3Dtiling(alphaBC,uBAR,NODESfaces,BasisUrb,ndim,COORref,DATAONLINE,DATALIM) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of RVE for each domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDOM = size(alphaBC,1) ;
[COLUMNS_RVEloc,nTYPErve,COLUMNS_RVE] = SeparatedModesOperations(DATAINM,nDOM) ;


% 
% if  
%     
% 
% DATAINM= DefaultField(DATAINM,'OBJECTIVE_FUNCION_WITH_STIFFNESS_MATRIX',0) ;
% 
%   NODES = cell2mat(NODESfaces) ;
%  NODES = unique(NODES) ;
%  DOFf = small2large(NODES,ndim) ;  % Boundary DOFf
%  DOFi = 1:(size(BasisUrb,1));  % Interior DOFf
%  DOFi(DOFf) = [] ;
% 
% 
% if  DATAINM.OBJECTIVE_FUNCION_WITH_STIFFNESS_MATRIX == 1
%   %  K = RecoverStiffnessMatrixUnitaryProject(DATAONLINE) ;
%     %   u = K(DOFi,DOFi)\ones(length(DOFi),1) ;
%     Kffbar = K(DOFf,DOFf) - K(DOFf,DOFi)*(K(DOFi,DOFi)\K(DOFi,DOFf)) ;
%     [INTTTG DOFremove DOFremove2] =  intersect(DOFf,DOFS_reference);
%     DOFf_loc = 1:length(DOFf) ;
%     DOFf_loc(DOFremove) = [] ;
%     Kffbar =     Kffbar(DOFf_loc,DOFf_loc)  ;
%     Hbar = chol(Kffbar) ;  % If its no positive definite, then check what happens with DOFS_reference
%     HbarINV = inv(Hbar) ;
%      kmax = (max(max(Kffbar))) ; 
%     HbarAUG = sparse(size(K,1),size(K,2)) ; 
%     HbarAUG(DOFf_loc,DOFf_loc) = Hbar ; 
%     for i = 1:length(DOFS_reference) % Add term to diagonal
%         HbarAUG(DOFS_reference(i),DOFS_reference(i) ) = sqrt(kmax) ; 
%     end
%     Hbar = HbarAUG ; 
%     
%      kmaxINV = 1/(max(max(Kffbar))) ; 
%     HbarAUGinv = sparse(size(K,1),size(K,2)) ; 
%     HbarAUGinv(DOFf_loc,DOFf_loc) = HbarINV ; 
%     for i = 1:length(DOFS_reference) % Add term to diagonal
%         HbarINV(DOFS_reference(i),DOFS_reference(i) ) = sqrt(kmaxINV) ; 
%     end
%     HbarINV = HbarAUGinv  ; 
%     
%     
% else
%     HbarINV  = [] ;
%     Hbar = [] ; 
% end
%  

% For later purposes, we define   global basis matrices
% ----------------------------------------------------
 [BasisUrbGLO,BasisUdefGLO,BasisRdefGLO] = GlobalBasisMatrices(BasisUrb,nDOM,BasisUdef,BasisRdef,COLUMNS_RVE);

% STEP 4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Assembly system of equations   A x = b
% % ----------------------------------------
% -----------------------------------
% A) Matrix Pcomp (compatibility matrix ), and vectors bCOMPr% qRB = PcompRR\bCOMPr ;
%  and bCOMPd
% % ---------------------------------
%[PcompRR PcompRD PcompDD bCOMPr bCOMPd] =
%AssemblyPcompNEW(BasisUrb,BasisUdef,f1,f2,nDOM,uBAR,ALPHA_ENDS) ;

 
tic
disp('Assembly P')
[Pcomp, bCOMP, mCOMP,INDrig,INDdef DATAwmethod] = ...
    AssemblyPcompTILING2D3D(BasisUrb,BasisUdef,NODESfaces,uBAR,alphaBC,COLUMNS_RVE,COLUMNS_RVEloc,ndim,...
    THETAfaces,GAMMAfaces,DATAINM) ;
toc
disp(['...Done  ' ])
% 
% kmax =  max(max(KdomRED{1})) ; 
% DATAINM = DefaultField(DATAINM,'FACTOR_MULTIPLY_DISPLACEMENT_EQUATIONS',0) ; 
% if DATAINM.FACTOR_MULTIPLY_DISPLACEMENT_EQUATIONS == 0
% FACTOR_P = 1;  
% else
%     FACTOR_P =  kmax*DATAINM.FACTOR_MULTIPLY_DISPLACEMENT_EQUATIONS ; 
% 
% end
% Pcomp = Pcomp*FACTOR_P ; 
% bCOMP = bCOMP*FACTOR_P ; 
% mCOMP = mCOMP*FACTOR_P ; 


 


% B) Rigid body reactions
% % RB reactions
[fextDOMglo,rRBglo,reactDOMrbGLO] = RigidBodyReactions(BasisUrb,nDOM,BasisRrb,FORCE_PROJECTS,INTENSITY_LOADS) ; 

%%%%%%%%%%%%%%%%%
% C) Transference matrix
% ----------------------
disp('Assembly of Treac')
tic
[Treac, cREAC, s, INDrigR,INDdefR]=  AssemblyTreacTILING2D3D(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
    THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb) ;
toc
disp('...Done')

DATAINM = DefaultField(DATAINM,'ADD_CONSTRAINT_RESULTANT_FORCES',0) ; 
if DATAINM.ADD_CONSTRAINT_RESULTANT_FORCES == 1
 
     [Tresult, cRESULT, sRESULT]=  AssemblyTreacRESULTANT(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
    THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb,COORref) ;
    
else
    Tresult = [] ; 
    cRESULT = [] ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D) Reaction-displacement basis matrix
% -----------------------------
 
[Hqr KdomREDglo] = RedStiffnessMatrix_Hmatrix(BasisUdef...
    ,nDOM,COLUMNS_RVE,KdomRED,BasisRdef,DATAINM,BasisUrb,BasisRrb) ;
 
%
% F) External forces
[Fbar,fextDOMglo,reactDOMrbGLO] = ...
    ExternalForcesReduced(reactDOMrbGLO,fextDOMglo,BasisUdefGLO,DATAINM,BasisUrbGLO,INDrig,INDdef) ; 

% G) Matrices for alternative formulations 
DATAwmethod = AlternativeFormulationsMatrices(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
    THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb,BasisUdef,BasisUrb,DATAwmethod,DATAINM,COORref,...
    INDrigR,INDdefR,INDrig,INDdef,uBAR) ; 


%% EQUATION
nRB = size(BasisUrb,2) ; 
DATAwmethod.INDrigR = INDrigR ;  DATAwmethod.INDdefR = INDdefR ;
disp('Solving linear system')
tic
[qRB,qDEF,rDEF] = DomainEquationAssembly(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
    KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAINM,DATAwmethod,rRBglo);
toc
disp('Done')
% %%%%%%%%%%%%%%%%%%%%%%

% Displacement field
d = BasisUrbGLO*qRB + BasisUdefGLO*qDEF ;
disp('MAx  y disp =')
disp(['MAxDISPY =',num2str(max(abs(d))'), ' (m)'])
REACT = reactDOMrbGLO + BasisRdefGLO*rDEF ;

% OBJECTIVE FUNCTION
%-------------------
% Contribution displacements
%
%dbstop('224')
[Q_d,Q_r] = ObjectiveFunction(qRB,Pcomp,INDrig,INDdef,qDEF,bCOMP,mCOMP,Treac,rDEF,cREAC,s) ; 

% % %See if global equilibrium  is fulfilled
% REACT_left_end = REACT(f1);
% resultantREAC = BasisUrb(f1,:)'*REACT_left_end;
% disp(['Reactions X = ',num2str(resultantREAC(1))]) ;
% disp(['Reactions Y = ',num2str(resultantREAC(2))]) ;
% if size(BasisUrb,2)==6
%     disp(['Reactions Z = ',num2str(resultantREAC(3))]) ;
%
% end
% Stress reconstruction
if exist('setPoints')==0
    setPoints = [];
end
DATAINM = DefaultField(DATAINM,'PRINT_AVERAGE_STRESSES_ON_ELEMENTS',1) ;
[stressGLO,setElementsRED] = ...
    StressReconstructionMULTI(qDEF,posgp,setPoints,BdomRED,Cglo,nDOM,DATAINM,Wdom,ndim,COLUMNS_RVE,COLUMNS_RVEloc) ;
DATA.setElementsRED = setElementsRED ;




%%%% GENERATING FINITE ELEMENT MESH by repeating unit cells
%[COORrecons,CNrecons,Materials] = MeshByRepetitionDOMnew(COORref,CNref,f1NOD,f2NOD,nDOM,MaterialType) ;
DOMAINS_TO_INCLUDE = [] ; 
[COORrecons,CNrecons,Materials] = MeshByRepetitionTILING2D3D(COORref,CNref,NODESfaces,DATAONLINE...
    ,MaterialType,DATAINM,nMAT,DOMAINS_TO_INCLUDE) ;


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
        %d = BasisUrbGLO*qRB + BasisUdefGLO*qDEF ;
    case 'RIGID'
        d = BasisUrbGLO*qRB;
    case 'STRAIN'
        d = BasisUdefGLO*qDEF ;
    otherwise
        error('Option not implemented')
end

strainGLO = [] ;
 DATA.MakeMeshByRepetition.nDOM = DATAONLINE.NdomX ; 
 
 DATA = DefaultField(DATA,'PostProcessWithNOSLICES',1) ; 
DATA = DefaultField(DATA,'MATERIAL',[]) ; 
 
 if DATA.PostProcessWithNOSLICES == 1 & ~isempty(DATA.MATERIAL)
    nmat = length(DATA.MATERIAL.PLY) ; 
    ndom = DATAONLINE.NdomX ; 
    TOTnmat = nmat*ndom ; 
    MAT_TypeORIG = (1:TOTnmat) ;
    MAT_TypeNOSLICES = repmat(((1:nmat)'),ndom,1) ; 
    NewMaterialType = zeros(size(Materials)); 
   
    for imat = 1:length(MAT_TypeORIG)
        INDLOC =find(Materials==imat) ; 
        NewMaterialType(INDLOC) = MAT_TypeNOSLICES(imat) ; 
    end
    
   Materials =  NewMaterialType ; 
end
 
 
GidPostProcess(COORrecons,CNrecons,TypeElement,d,strainGLO, ...
    stressGLO,  REACT,NAME_BASE_GIDfiles,posgp,NameFileMeshss,Materials,DATA);
%end
if ~isempty(setElementsRED)
    
    clipboard('copy',num2str(setElementsRED'))
    
end


dNORM = reshape(d,ndim,[]) ;
dNORM = sqrt(sum(dNORM.*dNORM,1)) ;
maxD = max(dNORM) ;


% QrGLO(iproj) = Q_r ;
% QdGLO(iproj)  = Q_d;
% maxDISP(iproj)  = max(abs(d)) ;




%
% figure(15)
% hold on
% xlabel('Modes')
% ylabel('MAX DISP (m)')
% h = plot(nDEFglo,maxDISP);
%
% figure(16)
% hold on
% xlabel('Modes')
% ylabel('Objective fun. Displacements')
% h = plot(nDEFglo,QdGLO);
%
% figure(17)
% hold on
% xlabel('Modes')
% ylabel('Objective fun. Reactions')
% h = plot(nDEFglo,QrGLO);
%
%
%
