function MinimizationDistanceReactionsMethod(NAMEWS,DATAONLINE,DATAINM)

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
    nMAT = length(MATERIAL.PLY) ;
    DATA.INPUTDATAfile = PROJECT_LOADS{iproject} ;
    DATA.NOCALCULATE_DISPLACEMENTS = 1 ;
    % Calling Finite Element elastostatic program (but only for computing external forces)
    DATAOUT = FE_ELASTOSTATIC(FUNinput,DATA) ;
    load(DATAOUT.nameWORKSPACE,'Fb','Ftrac','CN','COOR','MaterialType') ;
    
    %  load(DATAOUT.nameWORKSPACE,'Fb','Ftrac','CN','COOR','MaterialType','K') ;
    
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
                %                 disp('borrar esto')
                %           KdomRED{itype} =  BasisUdef{itype}'*(K*BasisUdef{itype}) ;
            end
        end
        
    end
    FORCE_PROJECTS{iproject} = Fb(IDXdofs) + Ftrac(IDXdofs) ;
end
% -------------------------------------------------------------------------
% STEP 3  - Determination interface nodes and Dirichlet Boundary conditions
% -------------------------------------------------------------------------
ndim  = size(COORref,2) ;
%[NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COORref,CNref,DATA) ;

% Determination of nodes pertaining to faces, edges and corners
[NODESbound] =  PointPlanesRBODY_GEN(COORref,CNref,DATA) ;

NODESfaces = NODESbound.PLANE ;

%NNFF = cellfun(@transpose,NODESfaces,'un',0);  % Transpose of all entries of the cell
DATAIN.RecalculateRigidBodyMatrix = 0;
if DATAIN.RecalculateRigidBodyMatrix  == 1
    f = unique(cell2mat(NODESfaces)) ;
    f = small2large(f,ndim);
    % Rigid body-reactions mode matrix
    BasisRrb = zeros(size(BasisUrb)) ;
    BasisRrb(f,:) = BasisUrb(f,:) ;
else
    % Taking all DOFs by default as restricted surfaces may give rise to
    % problems. In this case, it is better to retrieve it from memory
    load(NAMEWS,'BasisRrb');
end



% Dirichlet conditions
uBAR = DomainRedDirichledBC_2D3Dtiling(alphaBC,uBAR,NODESfaces,BasisUrb,ndim,COORref) ;

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
[Pcomp, bCOMP, mCOMP,INDrig,INDdef] = ...
    AssemblyPcompTILING2D3D(BasisUrb,BasisUdef,NODESfaces,uBAR,alphaBC,COLUMNS_RVE,COLUMNS_RVEloc,ndim,...
    THETAfaces,GAMMAfaces) ;

% B) Rigid body reactions
% % RB reactions
% % ---------------
nRB = size(BasisUrb,2) ; % number of rigid body modes
reactDOMrbGLO = cell(1,nDOM) ; % Cell containing the reactDOMrb variable
fextDOMglo = cell(1,nDOM) ;  % External forces applied on each subdomain
CovMixInv = inv(BasisUrb'*BasisRrb) ;
rRBglo = [] ;
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
    rRBglo = [rRBglo;rRB] ;
end
%%%%%%%%%%%%%%%%%
% C) Transference matrix
% ----------------------
disp('Assembly of Treac')
[Treac cREAC s]=  AssemblyTreacTILING2D3D(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
    THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb) ;

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
[qRB,qDEF,rDEF] = DomainEquationAssembly(nDOM,nRB,nDEF,nREAC,Pcomp, bCOMP,INDrig,INDdef,...
    KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAINM);
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
Q_d = (qRB'*Pcomp(INDrig,INDrig)*qRB + 2*qRB'*Pcomp(INDrig,INDdef)*qDEF+qDEF'*Pcomp(INDdef,INDdef)*qDEF...
    -2*qRB'*bCOMP(INDrig)-2*qDEF'*bCOMP(INDdef) + ...
    mCOMP) ;
Q_r = (rDEF'*Treac*rDEF +2*rDEF'*cREAC + s) ;
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