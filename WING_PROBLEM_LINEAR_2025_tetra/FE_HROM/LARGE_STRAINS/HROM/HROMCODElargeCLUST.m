function  [BasisUall_cl,RECONS_PK2stress_cl,OTHER_output ]=  HROMCODElargeCLUST(DATAFILE,PARAMETERS)
% HROM   CODE, large strains, local basis (clustering)
% JAHO, 21-Feb-2021  (adaptation of HROMCODElargeDYN)
if nargin ==0
    load('tmp.mat')
end
TIME_COMPUT=[] ;
format long g

% READING AND CONSTRUCTING INPUT DATA FOR THE FE CODE
% -----------------------------------------------------------
% Prototypical file, see
%/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/HROM/BEAM_2inpvar_SMALL_DISP/INPUTS_HROM_2param.m
disp('*************************************')
disp('Reading and constructing HROM input data (clustering)')
disp('*************************************')
TIME_COMPUT.INPUTS = tic ;
[BasisUall_cl,BasisStwo_cl,ECMdata_cl,...
    MESH,MATPRO_cl,OPERHROM_cl,Fbody_cl,Ftrac_cl,DISP_CONDITIONS_cl,INICOND,DATAHROM_cl,TransClustDATA,OTHER_output,DATAoffline]  = ...
    feval(DATAFILE,PARAMETERS) ;
TIME_COMPUT.INPUTS = toc(TIME_COMPUT.INPUTS) ;
disp(['...done in (',num2str(TIME_COMPUT.INPUTS ),' s)'])
disp('*************************************')






disp('*************************************')
disp('INITIALIZATIONS')
disp('*************************************')
% Variables at time n (and previous time steps if required. Initialization of
% snapshots )
OTHER_output.DISP_CONDITIONS_cl =DISP_CONDITIONS_cl ;
[DATAHROM_cl,VAR,SNAP] = INITIALIZATIONvar_cluster(DATAHROM_cl,INICOND)   ;%
%d = VAR_n.DISP(:,istep-1) ;
TIME_COMPUT.TIME_STEP_LOOP = tic ;
%DATAoffline.ILOADSTATES = PARAMETERS.ILOADSTATES;

[DATAHROM,CONVERGED,CLUSTER_SEQUENCE,ERROR_TRANSITION]=NONLINEARstaticCLUST(DATAHROM_cl,DISP_CONDITIONS_cl,VAR,OPERHROM_cl,SNAP,...
    Fbody_cl,Ftrac_cl,MATPRO_cl,TransClustDATA,DATAoffline) ;

OTHER_output.CLUSTER_SEQUENCE  =CLUSTER_SEQUENCE ;
OTHER_output.ERROR_TRANSITION  =ERROR_TRANSITION ;

TIME_COMPUT.TIME_STEP_LOOP = toc(TIME_COMPUT.TIME_STEP_LOOP) ;
OTHER_output.TIME_COMPUT = TIME_COMPUT;
disp(['...done in (',num2str(TIME_COMPUT.TIME_STEP_LOOP ),' s)'])

% ------------------------------------
% POST-PROCESS WITH GID
% RECONSTRUCTION DISPLACEMENTS
% -----------------------------
DATAinpGID = [];
DATAinpGID.OPERreconstr.DISP.BASIS =  BasisUall_cl ;
DATAinpGID.OPERreconstr.DISP.coeff = 1 ;
% RECONSTRUCTION PK2 STRESSES
% -----------------------------
DATAinpGID.OPERreconstr.PK2STRESS.BASIS =  BasisStwo_cl ;
DATAinpGID.OPERreconstr.PK2STRESS.coeff = cell(size(ECMdata_cl)) ;
for iclusters = 1:length(ECMdata_cl)
    DATAinpGID.OPERreconstr.PK2STRESS.coeff{iclusters} = ECMdata_cl{iclusters}.coeff_reconstr_STWO ;
end

RECONS_PK2stress_cl  =DATAinpGID.OPERreconstr.PK2STRESS ;


PARAMETERS = DefaultField(PARAMETERS,'PRINT_RIGID_BODY',1 ) ;
if PARAMETERS.PRINT_RIGID_BODY == 1
    DATAinpGID.ADDITION_VARIABLE.DISP = DISP_CONDITIONS_cl{1}.RIGID_BODY_MOTION ;
end


for icluster_STORE = 1:length(DATAHROM_cl.COMMON.STORE.NSTEPS_CLUSTER)
    GidPostProcessLARGE_cluster(MESH,DATAHROM_cl.COMMON,icluster_STORE,DATAinpGID,CLUSTER_SEQUENCE);
end






% External nodal forces
%
%     % Prescribed displacements
%
%
%
%     % ------------------------------------
%     % EXTERNAL ACTIONS AT THIS TIME STEP
%     % ------------------------------------
%     % A). External forces
%     % --------------------------
%     F_f = zeros(size(FORCES{1}.VALUE)) ;
%     for iforces = 1:length(FORCES) % Loop over set of forces (tractions, body ....)
%         TIME_FACTOR = FORCES{iforces}.TIME_FACTOR(istep);
%         F_f = F_f + FORCES{iforces}.VALUE*TIME_FACTOR;
%     end
%     % --------------------------------------------------
%     % B) Strain produced by Dirichlet Boundary conditions
%     % ---------------------------------------------
%     B_s_g0 = PRES_DISP.Bst_u*PRES_DISP.TIME_FACTOR(istep) ; %  BBnw_s*gBOUND(:,istep)  ;
%     %
%     increT = DATA.TIME_DISCRETIZATION(istep)- DATA.TIME_DISCRETIZATION(istep-1);
%     % ---------------------------
%     %  C) NEWTON-RAPHSON algorithm
%     %- ------------------------
%     %   tic
%     % NEWTON-RAPHSON
%     norm_res = 1;
%     tol_rel  = tolNEWTONRAPSHON;
%     iter = 1 ;
%     SALIR = 0 ;
%
%
%     DynDATA = DefaultField(DynDATA,'betaNM',1/4) ;
%     DynDATA = DefaultField(DynDATA,'gammaNM',1/2) ;
%
%
%     betaNM  = DynDATA.betaNM ;
%     gammaNM = DynDATA.gammaNM;
%
%     % Initial displacements (for iteration purposes)
%     Unp1 = Un ;
%     if  DynDATA.DYNAMIC == 0
%         U_k = Un(DOFf) ;
%     else
%         uTIL  = Un(DOFf) + VELn*increT +  0.5*(1-2*betaNM)*ACELn*increT^2 ;
%         U_k    =  uTIL ; %
%     end
%     deltaU = zeros(size(U_k)) ;
%
%     % if istep ==45
%     %     dbstop('34')
%     %     gammaNM
%     % end
%     Celas = [] ; % Elasticity matrix
%
%     TIMEelap.UpdateStress = 0 ;
%     TIMEelap.Residual = 0 ;
%     TIMEelap.AssesmblyK = 0 ;
%     TIMEelap.Solver = 0 ;
%     TIMEelap.Total = 0 ;
%
%     TimeTotal = tic;
%     while ((SALIR == 0) && (iter <= max_iter))
%
%
%
%         % dbstop('39')
%         eNP1 = BBfNW*U_k + B_s_g0  ; % Computing "strain-like" global vector field
%         % Computing vector of stresses and tangent moduli at each gauss point
%         % ------------------------------------------------------------------------
%
%         % dbstop('39')
%         CALC_CTANG = 1 ;
%
%         TIMEloc = tic;
%         if isempty(BBfT)
%             WEIGHTS = wST ;
%         else
%             WEIGHTS = [] ;
%         end
%
%
%         [Ctang,   stressSTNP1, EPnp1,  alphaNP1,  sigmayNP1,CtangM] = ...
%             UpdateStresses_J2_FE2D(eNP1,PROPMAT,EPn,sigmayN,alphaN,CALC_CTANG,iter,WEIGHTS) ;
%         TIMEloc = toc(TIMEloc) ;
%         TIMEelap.UpdateStress = TIMEelap.UpdateStress + TIMEloc ;
%
%
%
%         %     if STORE_CBMAT_SNAPSHOTS == 2
%         %         ifin = length(U_k)*iter; iini = length(U_k)*(iter-1)+1 ;
%         %         CBMAT(:,iini:ifin) = CtangNP1_B ;
%         %     end
%
%         % dbstop('52')
%
%         TIMEloc = tic;
%         % ----------------
%         if ~isempty(BBfT)
%             FINT = BBfT*stressSTNP1 ;
%         else
%             FINT = BBfNW'*(stressSTNP1.*wST) ;
%         end
%         TIMEloc = toc(TIMEloc) ;
%
%
%
%         %%%%
%         % External forces
%         if DynDATA.DYNAMIC == 0
%             FEXT = F_f ;
%             ACELnp1 = ACELn ;
%         else
%             % Dynamic case
%             % ------------
%             ACELnp1 =  (U_k-uTIL)/betaNM/increT^2;
%             FEXT = F_f -massMs*gDD- massMf*ACELnp1 ;
%         end
%         %%%%
%         %     if ~isempty(DATA_INTERPOLATORY)
%         %         if DATA_INTERPOLATORY.BLOCKDECOMP_EXPMAT == 1
%         %             FEXT = DATA_INTERPOLATORY.Dfext*FEXT ;
%         %         end
%         %     end
%         % Residual
%         ResF= FINT -FEXT;%  F_f;      % Computing residual vector
%         criter_f =  CheckConvergence(FINT,FEXT,ResF,num2str(iter));
%
%         TIMEelap.Residual = TIMEelap.Residual + TIMEloc ;
%
%
%
%         % -------------------------------------
%         if criter_f<tol_rel  | norm_res <tol_rel
%             SALIR = 1 ;
%             U_kp1 = U_k ;
%             %         % dbstop('78')
%             %         [ResF,FINT_Z] = CalF_RESmethod (ISRESID,BBfNWTnp,stressSTNP1,F_fNP,massMsNP,gDD,massMfNP,ACELnp1,...
%             %    DATA_FE1         nmodesU,CUBATURE_FINT_FEXT,FINT_Z,BBfNWT,FEXT,FINT) ;
%         else
%
%
%             TIMEloc = tic;
%
%             if isempty(ASSEMBLY_INFO)
%
%                 IMPLE = 0;
%
%                 if IMPLE  == 1
%                     % 25-June-2019 ---- Attempt to improve the efficiency ....
%                     Kstiff = BBfNW'*Ctang*BBfNW ;
%                     %             NPARTS = 100 ;
%                     %             INDgauss = unique(ceil(linspace(1,size(BBfNW,1),NPARTS))) ;
%                     %             Kstiff = sparse(size(BBfNW,2),size(BBfNW,2)) ;
%                     %
%                     %             for iparts = 1:NPARTS-1
%                     %                 IND = INDgauss(iparts):INDgauss(iparts+1) ;
%                     %                 Kstiff = Kstiff + BBfNW(IND,:)'*Ctang(IND,IND)*BBfNW(IND,:) ;
%                     %             end
%                 else
%
%                     % dbstop('92')
%                     if ~isempty(BBfT)
%                         CtangNP1_B = Ctang*BBfNW ;
%                         Kstiff = BBfT*CtangNP1_B ; % Jacobian matrix
%                     else
%                         Kstiff = Ctang*BBfNW ;
%                         Kstiff = BBfNW'*Kstiff ;     % Jacobian matrix  (Weights included in CtangNP1_B)
%                     end
%                 end
%             else
%                 Kstiff = AssemblyKMethodStandard_nonl(ASSEMBLY_INFO,CtangM) ;
%             end
%
%
%
%             if DynDATA.DYNAMIC ==1
%                 Kstiff = Kstiff + massMf/increT^2/betaNM  ;
%             end
%
%             TIMEloc = toc(TIMEloc) ;
%             TIMEelap.AssesmblyK = TIMEelap.AssesmblyK  + TIMEloc;
%
%             TIMEloc =tic ;
%             dirSEARCH = - Kstiff\ResF ;
%             TIMEloc = toc(TIMEloc) ;
%             TIMEelap.Solver = TIMEelap.Solver + TIMEloc ;
%
%             % Energ = -dirSEARCH'*ResF;
%
%             % if ~isempty(DATA_INTERPOLATORY) &&  ~isempty(DATA_INTERPOLATORY.ALPHA_PERTURB)
%
%             %             if Energ <0
%             %                 factorr = norm(Kstiff);
%             %                 alpha_per = DATA_INTERPOLATORY.ALPHA_PERTURB ;
%             %                 Kstiff = (alpha_per)*factorr*eye(size(Kstiff)) + (1-alpha_per)*Kstiff ;
%             %                 dirSEARCH = - Kstiff\ResF ;
%             %                 Energ = -dirSEARCH'*ResF;
%             %
%             %             end
%             %
%
%
%             %     end
%
%
%             U_kp1 = U_k + dirSEARCH;  % Displacement update
%             %dbstop('136')
%             %   deltaU =  U_kp1-U_k ;
%             iter = iter+1 ;
%             U_k = U_kp1 ;
%
%
%         end
%     end
%
%
%
%     % TESTSOLUTION =1 ;
%     % if TESTSOLUTION ==1
%     %
%     %
%     %     CtangNP1_B = Ctang*BBfNW ;
%     %              dbstop('143')
%     %             Kstiff = BBfT*CtangNP1_B ;
%     %
%     %     CBs = BBfT*Ctang*B_s_g0
%     %     dU = -Kstiff\CBs ;
%     % end
%
%
%
%
%
%     %%% Updating displacements and velocities
%     Unp1(DOFf) = U_kp1 ;
%
%     if DynDATA.DYNAMIC == 1
%         VEL_TILDE = VELn  + (1-gammaNM)*increT*ACELn;
%         VELnp1 =  VEL_TILDE + gammaNM*increT*ACELnp1;
%     else
%
%         VELnp1 = VELn ;
%     end
%
%     TimeTotal = toc(TimeTotal);
%
%     TIMEelap.Total = TimeTotal ;
%
%
%     %  if istep == 965
%     %         dbstop('158')FECODElargeDYN
%     %         Kstiff;
%     %
%     %
%     % end
%
%
%     % %%%
%     %   DETERMINE_EIGENFREQ = 1 ;
%     % % %
%     % if DETERMINE_EIGENFREQ == 1
%     %   %  if istep == 8000
%     %         dbstop('173')
%     %         disp('Eigenvalues')
%     %         K =   (BBfT*CtangNP1_B )  ;
%     %        % M =  (massMf) ;
%     %         [DIAGONALMval] = eigs(massMf,K);
%     %         DIAGONALMval = 1./DIAGONALMval ;
%     %         EIGENFREQ = sqrt((DIAGONALMval)) ;
%     %         EIGENPERIOD = 2*pi./EIGENFREQ ;
%     %         EIGENPERIOD = sort(EIGENPERIOD,'descend') ;
%     %         EIGENPERIOD
%     %
%     %
%     %    % end
%     % end
%
%     %dbstop('186')
%     % if STORE_CBMAT_SNAPSHOTS == 1
%     %     CBMAT = CtangNP1_B ;
%     % elseif STORE_CBMAT_SNAPSHOTS ==3
%     %     CBMAT = CtangNP1_B*deltaU ;
%     % end
%
%     CONVERGENCE =1 ;
%
%     if iter> max_iter
%         % dbstop('202')
%         warning('convergence error')
%         CONVERGENCE = 0 ;
%
%     end
%     % toc
%     %- Compute reactions
%     % ------------------
%     [NODESV_np1,DATA ]= ReactionsFE(BstW,GAUSSV_np1,FORCES,OPERfe,istep,NODESV_np1,wST,Bst,DATA) ;
%     % Compute Von-Mises stresses, element-wise
%     % ----------------------------------------
%     % See function StressStrains.m (to be implemented)
%     % -----------------------------------------
%     % Kinematic variables at slave DOFs
%     %.............................---------------
%     NODESV_np1 = UpdateKinematicVariables(OPERfe,PRES_DISP,BND_Dyn,NODESV_np1,istep) ;
%     % --------------------------
%     if CONVERGENCE == 0
%         istep_conv = istepSTORE ; % istep-1 ; % Last converged step
%         %         INDCONVER=  find(DATA.STEPS_TO_STORE == istep_conv ) ;
%         %
%         %         if   ~isempty(INDCONVER)
%         %             istep_conv = istepSTORE ;  % Last stored snapshot
%         %         end
%
%         %
%         %             istepSTORE = istepSTORE + 1;
%         %         end
%         %   break ;
%         %  pause(3)
%     end
%     % ----------------------------------------------------------
%     % VARIABLES UPDATE   (time n --> time np1)
%     % -----------------
%     [NODESV_n,GAUSSV_n] = UpdateHistoricalVariables(NODESV_np1,GAUSSV_np1) ;
%     % stressST, strain and displacement snapshot matrices
%     if    CONVERGENCE == 1
%
%         % REACTION FORCES ARE ALWAYS STORED IN MEMORY, sparse format (change 22-Jan-2020)
%
%         if find(DATA.STEPS_TO_STORE == istep )
%             istepSTORE = istepSTORE + 1;
%             [GAUSS_SNAP,NODES_SNAP] = UpdateSnapshots(GAUSSV_np1,NODESV_np1,istepSTORE,...
%                 GAUSS_SNAP,NODES_SNAP);
%         else
%
%         end
%         niter_TOT(istep) = iter ;
%     else
%         % Non-converged.
%         [GAUSS_SNAP,NODES_SNAP,DATA] = OnlyConvergedSteps(GAUSS_SNAP,NODES_SNAP,istep_conv,DATA);
%         break
%     end






