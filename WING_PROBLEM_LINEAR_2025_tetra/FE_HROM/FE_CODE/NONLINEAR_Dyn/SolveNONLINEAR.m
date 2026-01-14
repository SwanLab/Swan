function [NODES_SNAP,GAUSS_SNAP,OPERfe,NODESV_PROP,GAUSSV_PROP,DATA] = ...
    SolveNONLINEAR(Fb,Ftrac,dR,DOFr,COOR,CN,TypeElement,PROPMAT,...
    typePROBLEM,Bst,DOFm,Gb,DATA,wST,Nst,BstW,M,ndim,ASSEMBLY_INFO) ;

if nargin == 0
    load('tmp1.mat') 
end

% Copy of SolveElastic. Adaption to nonlinear problems
% -----------------------------------------------------
% This function returns   the (nnode*ndim x 1) vector of nodal displacements (d),
% as well as the arrays of stresses and strains
%%% points  (qheatGLO)
% K = Global stiffness matrix   (nnode*ndim x nnode*ndim)
% Fb = External force vector due to  body forces  (nnode*ndim x 1)
% Ftrac = External force vector due to  boundary tractions    (nnode*ndim x 1)
% DOFr = Set of restricted DOFs
% dR = Vector of prescribed displacements
% ----------------------
%dbstop('14')
if nargin == 0
    load('tmp0.mat')
end

nnode = size(COOR,1);  nelem = size(CN,1); nnodeE = size(CN,2) ;     %
[ngausT,DATA.ndof ]= size(Bst) ;
% -------------------------------------------------
% Degrees of freedom, number of elements ....
% -----------------------------------------------
OPERfe = FE_numbers_DOFS(nnode,ndim,DOFr,DOFm,Gb,nelem,DATA,nnodeE,ngausT) ;
OPERfe.ndimSP = size(COOR,2) ; 

% ......................................
% Operators BBf, BBfNW
[BBf,BBfNW,massMf,massMs] = BBf_BBfNW_Operators(OPERfe,Bst,BstW,M)  ;
% ----------------------------------------------
% ---------------------------------------------------------
% EXTERNAL ACTIONS (FORCES, AND PRESCRIBED DISPLACEMENTS )
% ------------------------------------------------------------
[FORCES,PRES_DISP,DATA] = ExtForces_Pdisp_TIME(DATA,dR,Fb,Ftrac,OPERfe,Bst) ;
% -------------------------------------------------------------
DynDATA.DYNAMIC = 0 ;  % For dynamic problems (so far dissabled)
% ------------------------------------------------------------
% Acceleration Dirich. Boundaries.
% -------------------------------------------------------------
BND_Dyn = AccelerationDirichletBoundary(DynDATA,OPERfe,DATA) ;
% ---------------------------------------------------------------

% Number of iterations
niter_TOT = zeros(OPERfe.nsteps,1);
% ------------------------------------------------------
% INITIALIZATION HISTORIC VARIABLES ("Nodes variables" and "Gauss variables", time n, and time np1)
% -----------------------------------------
DATA = DefaultField(DATA,'TypeImplementation','J2_plasticity_small_strains') ;
[NODESV_n,GAUSSV_n,NODESV_PROP,GAUSSV_PROP] = InitHistoricVariables(BND_Dyn,OPERfe,PROPMAT,DATA) ;
% --------------------
% SNAPSHOT MATRICES. Initialization
% --------------------
DATA = DefaultField(DATA,'STEPS_TO_STORE',1:OPERfe.nsteps) ; 
DATA.STEPS_TO_STORE = unique([1,DATA.STEPS_TO_STORE]) ; 

[NODES_SNAP,GAUSS_SNAP] = InitSnapShotMatrices(OPERfe,DynDATA,DATA,NODESV_n,GAUSSV_n,NODESV_PROP,GAUSSV_PROP)  ;
 
[DATA ] = ComputeReactionsONLINE(DATA,OPERfe)  ;% Created JAN-2020, to avoid storing reactions at all time steps 

% Loop over time steps
disp(['Loop over time steps...']) ;
Cglo = [] ; 
DATA = DefaultField(DATA,'STEPS_TO_STORE',1:length(DATA.TIME_DISCRETIZATION)) ; 
istepSTORE = 1; 
for istep = 2:length(DATA.TIME_DISCRETIZATION)
    disp(['Time step = ',num2str(istep)]) 
   
    % ------------------------------------
    % EXTERNAL ACTIONS AT THIS TIME STEP
    % ------------------------------------
    % A). External forces
    % --------------------------
    F_f = zeros(size(FORCES{1}.VALUE)) ;
    for iforces = 1:length(FORCES) % Loop over set of forces (tractions, body ....)
        TIME_FACTOR = FORCES{iforces}.TIME_FACTOR(istep);
        F_f = F_f + FORCES{iforces}.VALUE*TIME_FACTOR;
    end
    % --------------------------------------------------
    % B) Strain produced by Dirichlet Boundary conditions
    % ---------------------------------------------
    B_s_g0 = PRES_DISP.Bst_u*PRES_DISP.TIME_FACTOR(istep) ; %  BBnw_s*gBOUND(:,istep)  ;
    %
    increT = DATA.TIME_DISCRETIZATION(istep)- DATA.TIME_DISCRETIZATION(istep-1);
    % ---------------------------
    %  C) NEWTON-RAPHSON algorithm
    %- ------------------------
 %   tic
    [GAUSSV_np1,NODESV_np1,iter,CONVERGENCE] = ...
        NewtonRaphsonDyn_GEN(DATA,NODESV_n,GAUSSV_n,BBfNW,...
        istep,F_f,OPERfe,B_s_g0,...
        PROPMAT,massMf,massMs,BND_Dyn,DynDATA,increT,BBf,wST,ASSEMBLY_INFO) ;
  % toc
    %- Compute reactions 
    % ------------------
    [NODESV_np1,DATA ]= ReactionsFE(BstW,GAUSSV_np1,FORCES,OPERfe,istep,NODESV_np1,wST,Bst,DATA) ; 
    % Compute Von-Mises stresses, element-wise 
    % ----------------------------------------
    % See function StressStrains.m (to be implemented)
    % -----------------------------------------
    % Kinematic variables at slave DOFs
    %.............................---------------
    NODESV_np1 = UpdateKinematicVariables(OPERfe,PRES_DISP,BND_Dyn,NODESV_np1,istep) ;
    % --------------------------
    if CONVERGENCE == 0
        istep_conv = istepSTORE ; % istep-1 ; % Last converged step
%         INDCONVER=  find(DATA.STEPS_TO_STORE == istep_conv ) ; 
%         
%         if   ~isempty(INDCONVER)
%             istep_conv = istepSTORE ;  % Last stored snapshot 
%         end
        
%         
%             istepSTORE = istepSTORE + 1;
%         end
        %   break ;
        %  pause(3)
    end
    % ----------------------------------------------------------
    % VARIABLES UPDATE   (time n --> time np1)
    % -----------------
    [NODESV_n,GAUSSV_n] = UpdateHistoricalVariables(NODESV_np1,GAUSSV_np1) ;
    % stressST, strain and displacement snapshot matrices
    if    CONVERGENCE == 1
        
        % REACTION FORCES ARE ALWAYS STORED IN MEMORY, sparse format (change 22-Jan-2020)
        
        if find(DATA.STEPS_TO_STORE == istep )
            istepSTORE = istepSTORE + 1; 
            [GAUSS_SNAP,NODES_SNAP] = UpdateSnapshots(GAUSSV_np1,NODESV_np1,istepSTORE,...
                GAUSS_SNAP,NODES_SNAP);
        else
            
        end
        niter_TOT(istep) = iter ; 
    else
        % Non-converged. 
        [GAUSS_SNAP,NODES_SNAP,DATA] = OnlyConvergedSteps(GAUSS_SNAP,NODES_SNAP,istep_conv,DATA);
        break
    end
end

 

 