function [DATA_cl,CONVERGED,CLUSTER_SEQUENCE,ERROR_TRANSITION]=...
    NONLINEARstaticCLUST(DATA_cl,DISP_CONDITIONS_cl,VAR,OPERFE_cl,SNAP,Fbody_cl,Ftrac_cl,MATPRO_cl,...
    TransClustDATA,DATAoffline)

if nargin == 0
    load('tmp1.mat')
end

TransClustDATA = DefaultField(TransClustDATA,'SequenceClusters',[]) ;
% DETERMINATION  OF WHICH IS THE FIRST CLUSTER TO VISIT
% WE TAKE FOR GRANTED THAT THE INITIAL DISPLACEMENT IS ZERO
CLUSTER_SEQUENCE = zeros(1,length(DATA_cl.COMMON.STEPS)) ;
ERROR_TRANSITION = zeros(1,length(DATA_cl.COMMON.STEPS)) ;
DATAoffline = DefaultField(DATAoffline,'IdealIndexCluster',[]) ;


 
istep = 2;  icluster_STORE = 1;  % icluster_STORE IS FOR STORING PURPOSES


% SELECTION INITIAL CLUSTER
%**************************************************+
[DATAoffline,iCL] = ClusterInitialSelection...
    (DATAoffline,VAR,DISP_CONDITIONS_cl,Fbody_cl,Ftrac_cl,istep,TransClustDATA,DATA_cl,OPERFE_cl,MATPRO_cl) ;
% **************************************************

CLUSTER_SEQUENCE(1) = iCL ;

nsteps = length(DATA_cl.COMMON.STEPS) ;
disp('*************************************')
disp(['LOOP over time steps (NSTEPS =',num2str(nsteps),')'])
disp('*************************************')
%DISP_CONDITIONS = DISP_CONDITIONS_cl{iCL} ;
 DATA = DATA_cl.COMMON ; % Input data
 DATA = DefaultField(DATA,'SMALL_STRAIN_KINEMATICS',0) ; 
 DATA = DefaultField(DATA,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ;

  DATA = DefaultField(DATA,'SNAP_ITER',[]) ; 
  DATA.NEWTON_RAPHSON = DefaultField(DATA.NEWTON_RAPHSON,'ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS',0) ; 
  DATA.INTERNAL_FORCES_USING_precomputed_Kstiff = 0 ; 
  
  if ~isstruct(DATA.STORE.TOLERANCE_SVD_COMPRESSION)
      TOL = DATA.STORE.TOLERANCE_SVD_COMPRESSION ; 
      DATA.STORE.TOLERANCE_SVD_COMPRESSION = [] ; 
      fff = fieldnames(DATA.STORE.VAR) ; 
          
      for  iii = 1:length(fff)
          DATA.STORE.TOLERANCE_SVD_COMPRESSION.(fff{iii}) =  TOL ;
      end
  %    DATA.STORE.TOLERANCE_SVD_COMPRESSION.DISP = 0 ;
      
  end
  
  
  

while istep <=nsteps
    
%         if  istep == 25
%              disp('borrar esto...')
%          end
%     if   ~isempty(TransClustDATA.SequenceClusters)
%         iCL = TransClustDATA.SequenceClusters(istep) ;
%     end
    
    CLUSTER_SEQUENCE(istep) = iCL ;    % Current cluster
    % Degrees of freedom
    
%     if iCL == 553
%         disp('borrar esto')
%     end
   
    disp('--------------------------------------------------')
    disp(['istep = ',num2str(istep), '  t_i =',num2str(DATA.STEPS(istep)),'  CLUSTER = ',num2str(iCL)])
    disp('---------------------------------------------------')
    % 1.a) Initialization displacements
    % Prescribed macroscopic deformation at  time t_n+1
    VAR = BoundaryConditionsImpose(VAR,DISP_CONDITIONS_cl,Fbody_cl,Ftrac_cl,iCL,istep) ;
    % *********************************
    %%%% NEWTON-RAPHSON ALGORITHM
    % *****************************
    OPERFE_cl{iCL}.HYDRO = [] ; OPERFE_cl{iCL}.KinternalFORCES_given = [] ; OPERFE_cl{iCL}.DOFm = [] ; 
    [VAR,CONVERGED ]= NewtonRapshonStaticLarge(DATA,OPERFE_cl{iCL},VAR,MATPRO_cl{iCL},DISP_CONDITIONS_cl{iCL}.DOFl) ;
    % ***************
    % Store snaphots
    % ***************
    [SNAP,DATA,icluster_STORE] = StoreInfoSnapshots(istep,icluster_STORE,SNAP,VAR,DATA,CONVERGED) ;
    if CONVERGED == 0
        disp('Convergence error ....')
        break
    end
    % ******************************
    % Transition between clusters
    % ******************************
   % if   isempty(TransClustDATA.SequenceClusters)
   if size(TransClustDATA.TransMatrix,1) == size(TransClustDATA.TransMatrix,2)
       % Version before July 24th 2022 (Not efficient in terms of memory management )
        [iCL,VAR,ERROR_TRANSITION(istep)]= TransitionClusterHROM(iCL,TransClustDATA,VAR,DATAoffline,DATA,...
            OPERFE_cl,MATPRO_cl,DISP_CONDITIONS_cl,istep,Fbody_cl,Ftrac_cl) ;
   else
       % % Version after July 24th 2022 (More efficient in terms of memory management )
        [iCL,VAR,ERROR_TRANSITION(istep)]= TransitionClusterHROMeff(iCL,TransClustDATA,VAR,DATAoffline,DATA,...
            OPERFE_cl,MATPRO_cl,DISP_CONDITIONS_cl,istep,Fbody_cl,Ftrac_cl) ;
   end
        
   % end
    istep = istep + 1;
end

