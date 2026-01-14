clc
clear all
% Structural analysis using domain-wise reduced-order data
% Prototype for one single domain
% See Domain_Decom_SVD.m
% Online stage. For the offline stage see DomainPartitionModal.m
%---------------------------------------------------------------
% INPUTS
%-----------------------------------------------------------
if exist('INPUT_PERIODIC')==0
    addpath('FE_CODE') ;
end
addpath('DATA_input')
INPUT_DATAFILE  ='DATA_BEAM3Drep'; 'DATAStent1' ; 'DATA_BEAMhexaIND'; 'DATA_BEAMhexa';  'DATA_BEAMporousREP' ;'DATA_BEAMsolidEMP30' ;  'DATA_BEAMslice';

DATA.NAMEWS = ['DATAWS/',INPUT_DATAFILE,'_OFFLINE','.mat'] ;
outputFE = ['DATAWS/',INPUT_DATAFILE,'_WS.mat'] ;
nDEFglo = [5:10];
nDOM = 40; % Number of domains
DOMAINS_with_FORCE = zeros(nDOM,1) ;
DOMAINS_with_FORCE(:)= 1 ;


eval(INPUT_DATAFILE)
for iproj = 1:length(nDEFglo)
    nDEF = nDEFglo(iproj) ;
    
    %  Domains subjected to the distributed force  (boolean)
    % Number of deformation modes
    nREAC =nDEF ;
    %     DATA.TypeUnitCell = 'HEXAG_2D_SQUARE';
    %     DATA = Defaiu
    
    % END INPUTS
    % ----------------------------------------------------------
    nDEF = min(nDOM,nDEF);
    nREAC = min(nREAC,nDEF) ;
    
    
    % Retrieval of reference mesh, coordinates, modes...
    %---------------------------------------------------
    load(outputFE,'Bst','wSTs','Cglo','Fpnt','CNb','Tnod','TypeElementB','CONNECTb','K');
    load(DATA.NAMEWS,'U','S','V','CNref','NODESref','COOR','TypeElement','posgp','BasisRdef','ElementsREF')
    % Computation of basis matrices, reduced stiffness matrix and external
    % forces
    [COORdom,BasisUrb,BasisUdef,BasisRdef,BasisRrb,f1,f2,fextDOM,KdomRED,f1NOD,f2NOD] ...
        = BasisMatricesRefDomain(COOR,NODESref,CNref,U,DATA,nDEF,nREAC,BasisRdef,posgp,...
        ElementsREF,Cglo,Bst,Fpnt,CNb,Tnod,TypeElementB,CONNECTb) ;
    ndim = size(COOR,2) ;
    if ndim == 2
        uBAR_0 = [0 0 0]; % Translation x, Translation y, Rotation around bottom corner (left end)
    else
        uBAR_0 = [0 0 0 0 0 0 ];
    end
    % Dirichlet conditions left end
    % ------------------------------
    if ndim ==2
        uBAR = zeros(size(f1)) ;
        ndim = size(COOR,2) ;
        uBAR(1:2:end) = uBAR_0(1) ;
        uBAR(2:2:end) = uBAR_0(2) ;
        uBAR= uBAR + uBAR_0(3)*BasisUrb(f1,3) ;
    else
        uBAR = zeros(size(f1)) ;
        
        uBAR(1:3:end) = uBAR_0(1) ;
        uBAR(2:3:end) = uBAR_0(2) ;
        uBAR(3:3:end) = uBAR_0(3) ;
        uBAR= uBAR + uBAR_0(4)*BasisUrb(f1,4) ;
        uBAR= uBAR + uBAR_0(5)*BasisUrb(f1,5) ;
        uBAR= uBAR + uBAR_0(6)*BasisUrb(f1,6) ;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Assembly system of equations   A x = b
    % ----------------------------------------
    
    %%
    % Matrix Pcomp (compatibility matrix ), and vectors bCOMPr% qRB = PcompRR\bCOMPr ;
    % and bCOMPd
    % ---------------------------------
    [PcompRR PcompRD PcompDD bCOMPr bCOMPd] =  AssemblyPcomp(BasisUrb,BasisUdef,f1,f2,nDOM,uBAR) ;
    
    %Checking that pure rigid body motions can be seamlessly reproduced
    % qRB = PcompRR\bCOMPr ;
    
    % For later purposes, we define the global basis matrices of displaceemnts
    BasisUrbGLO = repmat(sparse(BasisUrb),1,nDOM) ;
    BasisUrbGLO = mat2cell(BasisUrbGLO,size(BasisUrb,1),repmat(size(BasisUrb,2),nDOM,1)) ;
    BasisUrbGLO = blkdiag(BasisUrbGLO{:}) ;
    
    BasisUdefGLO = repmat(sparse(BasisUdef),1,nDOM) ;
    BasisUdefGLO = mat2cell(BasisUdefGLO,size(BasisUdef,1),repmat(size(BasisUdef,2),nDOM,1)) ;
    BasisUdefGLO = blkdiag(BasisUdefGLO{:}) ;
    
    BasisRdefGLO = repmat(sparse(BasisRdef),1,nDOM) ;
    BasisRdefGLO = mat2cell(BasisRdefGLO,size(BasisRdef,1),repmat(size(BasisRdef,2),nDOM,1)) ;
    BasisRdefGLO = blkdiag(BasisRdefGLO{:}) ;
    % dRB = BasisUrbGLO*qRB ;
    % ---------------------------------------------------------------
    %%
    % ---------------
    % RB reactions
    % ---------------
    nRB = size(BasisUrb,2) ;
    reactDOMrbGLO = cell(1,nDOM) ;
    fextDOMglo = cell(1,nDOM) ;
    %dbstop('107')
    CovMixInv = inv(BasisUrb'*BasisRrb) ;
    for idom = 1:nDOM
        fextDOMglo{idom} = DOMAINS_with_FORCE(idom)*fextDOM ;
        fextDOMrb = BasisUrb'*fextDOMglo{idom} ;  % Resultant (force x, force y and moment around reference point) of external forces
        rRB =- CovMixInv*fextDOMrb ;
        reactDOMrbGLO{idom} = BasisRrb*rRB ;
    end
    
    [Treac cREAC s]=  AssemblyTreac(BasisRdef,f1,f2,nDOM,reactDOMrbGLO) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reaction-displacement basis matrix
    % -----------------------------
    COV =BasisUdef'*BasisRdef ;
    Hqr= repmat(sparse(COV),1,nDOM) ;
    Hqr = mat2cell(Hqr,size(Hqr,1),repmat(size(COV,2),nDOM,1)) ;
    Hqr = blkdiag(Hqr{:}) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stiffness matrix
    % -----------------------------
    KdomREDglo= repmat(sparse(KdomRED),1,nDOM) ;
    KdomREDglo = mat2cell(KdomREDglo,size(KdomREDglo,1),repmat(size(KdomRED,2),nDOM,1)) ;
    KdomREDglo = blkdiag(KdomREDglo{:}) ;
    
    % External forces
    reactDOMrbGLO = cell2mat(reactDOMrbGLO) ;
    reactDOMrbGLO = (reactDOMrbGLO(:)) ;
    fextDOMglo  = cell2mat(fextDOMglo) ;
    fextDOMglo = (fextDOMglo(:)) ;
    
    Fbar = BasisUdefGLO'*(fextDOMglo+reactDOMrbGLO) ;
    
    %% EQUATION
    [A,b] = DomainEquationAssembly(nDOM,nRB,nDEF,nREAC,PcompRR,PcompRD,PcompDD,...
        KdomREDglo,bCOMPd,bCOMPr,Treac,Hqr,Fbar,cREAC);
    %%%%%%%%%%%%%%%%%%%%%%
    
    %%% SOLUTION
    x = full(A\b) ;
    %
    iacum = 0 ;
    qRB = x(1:nDOM*nRB) ;
    iacum =  nDOM*nRB ;
    qDEF = x((iacum+1):(iacum+nDOM*nDEF)) ;
    iacum = iacum+nDOM*nDEF ;
    rDEF = x((iacum+1):(iacum+nDOM*nREAC)) ;
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
    Q_d = qRB'*PcompRR*qRB + 2*qRB'*PcompRD*qDEF+qDEF'*PcompDD*qDEF-2*qRB'*bCOMPr-2*qDEF'*bCOMPd + norm(uBAR)^2 ;
    Q_r = rDEF'*Treac*rDEF +2*rDEF'*cREAC + s ;
    disp('--------------------------------------')
    disp(['Objective function DISPLACEMENT'])
    disp(['Q_q =',num2str((Q_d))])
    disp('--------------------------------------')
    disp(['Objective function REACTIONS'])
    disp(['Q_q =',num2str((Q_r))])
    
    % %See if global equilibrium  is fulfilled
    REACT_left_end = REACT(f1);
    resultantREAC = BasisUrb(f1,:)'*REACT_left_end;
    disp(['Reactions X = ',num2str(resultantREAC(1))]) ;
    disp(['Reactions Y = ',num2str(resultantREAC(2))]) ;
    if size(BasisUrb,2)==6
        disp(['Reactions Z = ',num2str(resultantREAC(3))]) ;
        
    end
    
    
    
    
    %%%% GENERATING FINITE ELEMENT MESH by repeating unit cells
    [COORrecons,CNrecons,Materials] = MeshByRepetitionDOM(COORdom,CNref,NODESref,f1NOD,f2NOD,nDOM) ;
    % -----------------------------------------------------------------
    
    
    
    
    
    % -----------------------------------------
    %%% PRINTING GID RESULTS
    % -----------------------
    if iproj == length(nDEFglo)
        NAME_BASE_GIDfiles = [INPUT_DATAFILE,'_APPROX','_',num2str(nDEF)] ;
        NameFileMeshss = [] ;
        strainGLOgid = [] ;
        stressGLOgid = [] ;
        GidPostProcess(COORrecons,CNrecons,TypeElement,d,strainGLOgid, ...
            stressGLOgid,  REACT,NAME_BASE_GIDfiles,posgp,NameFileMeshss,Materials,DATA);
    end
    
    QrGLO(iproj) = Q_r ;
    QdGLO(iproj)  = Q_d;
    maxDISP(iproj)  = max(abs(d)) ;
    
end


figure(15)
hold on
xlabel('Modes')
ylabel('MAX DISP (m)')
h = plot(nDEFglo,maxDISP);

figure(16)
hold on
xlabel('Modes')
ylabel('Objective fun. Displacements')
h = plot(nDEFglo,QdGLO);

figure(17)
hold on
xlabel('Modes')
ylabel('Objective fun. Reactions')
h = plot(nDEFglo,QrGLO);



