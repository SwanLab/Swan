function [Q_d,Q_r,maxD] =OnlineDomainSVD_MULTI(NAME_PROJECT_TEST,PROJECT_LOADS,INTENSITY_LOADS...
    ,NAMEWS,ALPHA_ENDS,uBAR_ENDS,DATAINM)
% Structural analysis using domain-wise reduced-order modeling
% Prototype for one single domain
% See MultiLearn.pdf
%---------------------------------------------------------------
% INPUTS
%-----------------------------------------------------------
%dbstop('10')
if nargin == 0
    load('tmp1.mat')
end
if exist('INPUT_PERIODIC')==0 ;     addpath('FE_CODE') ; end
if exist(PROJECT_LOADS{1}) ==0 ;     addpath('DATA_input'); end
% ----------
% STEP 1
% ---------
% Retrieving OFFLINE information
load(NAMEWS,'BasisUdef','BasisRdef','BasisUrb','CNref',...
    'COORref','NODESfacesLOC','TypeElement','posgp','Wdom','BdomRED') ;
DATAINM = DefaultField(DATAINM,'CUBATURE',[]) ;
DATAINM.CUBATURE = DefaultField(DATAINM.CUBATURE,'ACTIVE',0) ;

if DATAINM.CUBATURE.ACTIVE == 1
    load(NAMEWS,'BasisS','BdomRED','setPoints','WdomRED')
end
nDEF = size(BasisUdef{1},2) ; nREAC = size(BasisRdef{1},2) ;
% ----------
% STEP 2
% ---------
% LOOP OVER PROJECTS DEFINING EXTERNAL FORCES AND STIFFNESS MATRIX
% -----------------------------------------------------------------
%
FORCE_PROJECTS = {} ;
for iproject = 1:length(PROJECT_LOADS)
    eval(PROJECT_LOADS{iproject}) ;
    DATA.INPUTDATAfile = PROJECT_LOADS{iproject} ;
    DATA.NOCALCULATE_DISPLACEMENTS = 1 ;
    % Calling Finite Element elastostatic program (but only for computing external forces)
    DATAOUT = FE_ELASTOSTATIC(FUNinput,DATA) ;
    load(DATAOUT.nameWORKSPACE,'Fb','Ftrac','CN','COOR','MaterialType') ;
    
    [IDX D]= knnsearch(COOR,COORref) ;
    
    if any(abs(D) >1e-16 )
        dbstop('51')
        error('Non-conforming meshes')
    end
    IDXdofs = Nod2DOF(IDX,size(COOR,2)) ;
    % All properties are set in terms of the numbering of reference mesh
    %  dbstop('55')
    if iproject == 1
        load(DATAOUT.nameWORKSPACE,'Cglo') ;
        %        K_orig = K(IDXdofs,IDXdofs) ;
        %        KdomRED_orig = BasisUdef'*K*BasisUdef ;
        %%% NEW METHOD
        %        nstrain = size(BdomRED,1)/length(Wdom) ;
        %        wDIAG = CompWeightDiag(Wdom,nstrain)  ;
        %   dbstop('63')
        
        %
        KdomRED = {} ;
        for itype = 1:length(BdomRED)
            nstrain = size(BdomRED{1},1)/length(Wdom) ;
            
            if  DATAINM.CUBATURE.ACTIVE == 1
                setIndices =  small2large(setPoints{itype},nstrain) ;
                % Cglo is multiplied by Wdom(setPoints). Accordingly, we
                % define
                rW = WdomRED{itype}./Wdom(setPoints{itype}) ;
                % And make
                rwDIAG = CompWeightDiag(rW,nstrain)  ;
                KdomRED{itype} = (rwDIAG*BdomRED{itype}(setIndices,:))'*(Cglo(setIndices,setIndices)*BdomRED{itype}(setIndices,:)) ;
                
            else
                KdomRED{itype} =   BdomRED{itype}'*(Cglo*BdomRED{itype}) ;
            end
        end
        
    end
    FORCE_PROJECTS{iproject} = Fb(IDXdofs) + Ftrac(IDXdofs) ;
end
% -------------------------------------------------------------------------
% STEP 3  - Determination interface nodes and Dirichlet Boundary conditions
% -------------------------------------------------------------------------
ndim  = size(COORref,2) ;
[NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COORref,CNref,DATA) ;
f1 = small2large(f1NOD,ndim) ; % Degrees of freedom face F1 (local)
f2 = small2large(f2NOD,ndim) ; % Degrees of freedom face F2 (local)
f = [f1; f2] ;
BasisRrb = zeros(size(BasisUrb)) ;
BasisRrb(f,:) = BasisUrb(f,:) ;

% Dirichlet conditions
uBAR = DomainRedDirichletConditionss(uBAR_ENDS,ALPHA_ENDS,f1,f2,BasisUrb,ndim) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of RVE for each domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDOM = length(INTENSITY_LOADS{1}) ;
DATAINM = DefaultField(DATAINM,'SEPARATED_MODES',[]) ; % = TYPE_RVE ;
DATAINM.SEPARATED_MODES = DefaultField(DATAINM.SEPARATED_MODES,'COLUMNS_RVETYPE_TEST',[]) ; % = TYPE_RVE ;
if ~isempty(DATAINM.SEPARATED_MODES.COLUMNS_RVETYPE_TEST) && DATAINM.SEPARATED_MODES.ACTIVE ==1
    nTYPErve = length(DATAINM.SEPARATED_MODES.COLUMNS_RVETYPE_TEST) ;
    COLUMNS_RVEloc = {ones(1,nDOM)} ;
    COLUMNS_RVE = cell(1,nTYPErve) ;
    iacumPROJ = 0 ;
    [COLUMNS_RVE,COLUMNS_RVEloc ]=...
        ColumnsSeparatedRVEs(DATAINM,nTYPErve,nDOM,COLUMNS_RVE,iproject,iacumPROJ,COLUMNS_RVEloc);
    COLUMNS_RVEloc = COLUMNS_RVEloc{1} ;
else
    COLUMNS_RVEloc = ones(1,nDOM) ;
    nTYPErve = 1 ;
    COLUMNS_RVE = {1:nDOM} ;
end

% For later purposes, we define   global basis matrices
% ----------------------------------------------------
BasisUrbGLO = repmat(sparse(BasisUrb),1,nDOM) ;
BasisUrbGLO = mat2cell(BasisUrbGLO,size(BasisUrb,1),repmat(size(BasisUrb,2),nDOM,1)) ;
BasisUrbGLO = blkdiag(BasisUrbGLO{:}) ;

% Deformation basis matrices
BasisUdefGLO = DiagonalGlobalMatrixRVEs(nDOM,BasisUdef,COLUMNS_RVE) ;
% Reaction basis matrices
BasisRdefGLO = DiagonalGlobalMatrixRVEs(nDOM,BasisRdef,COLUMNS_RVE) ;
%
% BasisRdefGLO = repmat(sparse(BasisRdef),1,nDOM) ;
% BasisRdefGLO = mat2cell(BasisRdefGLO,size(BasisRdef,1),repmat(size(BasisRdef,2),nDOM,1)) ;
% BasisRdefGLO = blkdiag(BasisRdefGLO{:}) ;

% STEP 4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Assembly system of equations   A x = b
% % ----------------------------------------
% -----------------------------------
% A) Matrix Pcomp (compatibility matrix ), and vectors bCOMPr% qRB = PcompRR\bCOMPr ;
%  and bCOMPd
% % ---------------------------------
%[PcompRR PcompRD PcompDD bCOMPr bCOMPd] =  AssemblyPcompNEW(BasisUrb,BasisUdef,f1,f2,nDOM,uBAR,ALPHA_ENDS) ;
[PcompRR PcompRD PcompDD bCOMPr bCOMPd] = ...
    AssemblyPcompMULTI(BasisUrb,BasisUdef,f1,f2,nDOM,uBAR,ALPHA_ENDS,COLUMNS_RVE,COLUMNS_RVEloc) ;

% B) Rigid body reactions
% % RB reactions
% % ---------------
nRB = size(BasisUrb,2) ; % number of rigid body modes
reactDOMrbGLO = cell(1,nDOM) ; % Cell containing the reactDOMrb variable
fextDOMglo = cell(1,nDOM) ;  % External forces applied on each subdomain
CovMixInv = inv(BasisUrb'*BasisRrb) ;
for idom = 1:nDOM
    for iproj = 1:length(FORCE_PROJECTS)
        if isempty(fextDOMglo{idom})
            fextDOMglo{idom} = FORCE_PROJECTS{iproj}*INTENSITY_LOADS{iproj}(idom) ;
        else
            fextDOMglo{idom} =fextDOMglo{idom} +  FORCE_PROJECTS{iproj}*INTENSITY_LOADS{iproj}(idom) ;
        end
    end
    fextDOMrb = BasisUrb'*fextDOMglo{idom} ;  % Resultant (force x, force y and moment around reference point) of external forces
    rRB =- CovMixInv*fextDOMrb ;
    reactDOMrbGLO{idom} = BasisRrb*rRB ;
end
%%%%%%%%%%%%%%%%%
% C) Transference matrix
% ----------------------
%[Treac cREAC s]=  AssemblyTreacNEW(BasisRdef,f1,f2,nDOM,reactDOMrbGLO,ALPHA_ENDS) ;
[Treac cREAC s]=  AssemblyTreacMULTI(BasisRdef,f1,f2,nDOM,reactDOMrbGLO,ALPHA_ENDS,COLUMNS_RVE,COLUMNS_RVEloc) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D) Reaction-displacement basis matrix
% -----------------------------
% COV =BasisUdef'*BasisRdef ;
% Hqr= repmat(sparse(COV),1,nDOM) ;
% Hqr = mat2cell(Hqr,size(Hqr,1),repmat(size(COV,2),nDOM,1)) ;
% Hqr = blkdiag(Hqr{:}) ;

COV = cell(size(BasisUdef)) ;
for itype = 1:length(BasisUdef)
    COV{itype} =BasisUdef{itype}'*BasisRdef{itype} ;
end
Hqr = DiagonalGlobalMatrixRVEs(nDOM,COV,COLUMNS_RVE) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E) Stiffness matrix
% -----------------------------
KdomREDglo = DiagonalGlobalMatrixRVEs(nDOM,KdomRED,COLUMNS_RVE) ;

% KdomREDglo= repmat(sparse(KdomRED),1,nDOM) ;
% KdomREDglo = mat2cell(KdomREDglo,size(KdomREDglo,1),repmat(size(KdomRED,2),nDOM,1)) ;
% KdomREDglo = blkdiag(KdomREDglo{:}) ;
%
% F) External forces
reactDOMrbGLO = cell2mat(reactDOMrbGLO) ;
reactDOMrbGLO = (reactDOMrbGLO(:)) ;
fextDOMglo  = cell2mat(fextDOMglo) ;
fextDOMglo = (fextDOMglo(:)) ;

Fbar = BasisUdefGLO'*(fextDOMglo+reactDOMrbGLO) ;

%% EQUATION
[A,b] = DomainEquationAssembly(nDOM,nRB,nDEF,nREAC,PcompRR,PcompRD,PcompDD,...
    KdomREDglo,bCOMPd,bCOMPr,Treac,Hqr,Fbar,cREAC);
% %%%%%%%%%%%%%%%%%%%%%%

%%% SOLUTION
disp(['Solving reduced order system (number of eq. =',length(b),')'])
tic
x = full(A\b) ;
toc
%
iacum = 0 ;
qRB = x(1:size(PcompRR,1)) ;
iacum =  size(PcompRR,1);
qDEF = x((iacum+1):(iacum+ size(PcompDD,2))) ;
iacum = iacum+size(PcompDD,2);
rDEF = x((iacum+1):(iacum+size(Treac,2))) ;
%%%%%
% Displacement field
d = BasisUrbGLO*qRB + BasisUdefGLO*qDEF ;
disp('MAx  y disp =')
disp(['MAxDISPY =',num2str(max(abs(d))'), ' (m)'])
REACT = reactDOMrbGLO + BasisRdefGLO*rDEF ;

% OBJECTIVE FUNCTION
%-------------------
% Contribution displacements
%
Q_d = qRB'*PcompRR*qRB + 2*qRB'*PcompRD*qDEF+qDEF'*PcompDD*qDEF-2*qRB'*bCOMPr-2*qDEF'*bCOMPd + ...
    ALPHA_ENDS(1)*norm(uBAR{1})^2 +  ALPHA_ENDS(end)*norm(uBAR{end})^2 ;
Q_r = rDEF'*Treac*rDEF +2*rDEF'*cREAC + s ;
disp('--------------------------------------')
disp(['Objective function DISPLACEMENT'])
disp(['Q_q =',num2str((Q_d))])
disp('--------------------------------------')
disp(['Objective function REACTIONS'])
disp(['Q_q =',num2str((Q_r))])

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

[stressGLO,setElementsRED] = ...
    StressReconstructionMULTI(qDEF,posgp,setPoints,BdomRED,Cglo,nDOM,DATAINM,Wdom,ndim,COLUMNS_RVE,COLUMNS_RVEloc) ;
DATA.setElementsRED = setElementsRED ;




%%%% GENERATING FINITE ELEMENT MESH by repeating unit cells
[COORrecons,CNrecons,Materials] = MeshByRepetitionDOMnew(COORref,CNref,f1NOD,f2NOD,nDOM,MaterialType) ;
% -----------------------------------------------------------------





% -----------------------------------------
%%% PRINTING GID RESULTS
% -----------------------
%if iproj == length(nDEFglo)
NAME_BASE_GIDfiles = [NAME_PROJECT_TEST,'_APPROX','_',num2str(nDEF)] ;
NameFileMeshss = [] ;
ndim =  size(COORrecons,2) ;

strainGLO = [] ;
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
