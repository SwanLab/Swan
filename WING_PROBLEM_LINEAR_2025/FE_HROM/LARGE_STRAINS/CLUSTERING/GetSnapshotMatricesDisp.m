function [BasisU,SNAPdisp,MATRIX_POINTS_SPACE_PARAMETER,MESH,TIME_STEPS,TESTING_PARAM,INDEXES_TRAJECTORIES,DISP_CONDITIONS,...
    OPERFE,MATPRO,DATA]= GetSnapshotMatricesDisp(CASES,NAMEsnap_base,DATAoffline)


% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;

%  if    DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
%      % We are interested in a global basis matrix, also including the
%      % constrained DOFs
%      SNAPdisp_ALL =cell(1,length(CASES)) ;
%  end
TIME_STEPS = cell(size(CASES));

TimeDiscretizationLocal = TimeDiscretization ;


%DATAoffline.FunctionInputForces_original = 'InputForces' ; 
DATAoffline  =DefaultField(DATAoffline,'FunctionInputForces_original','InputForces') ; 
%DATAoffline  =DefaultField(DATAoffline,'NormalizedTrainingSnapshots',0) ; 

 
 INPUT_PARAMETERS = feval(DATAoffline.FunctionInputForces_original) ; %InputForces;

%DATAoffline  =DefaultField(DATAoffline,'IndexesInputForcesSubSampling',1:length(INPUT_PARAMETERS)) ; 

INPUT_PARAMETERS = INPUT_PARAMETERS(CASES) ; 

nstepsLOC = length( TimeDiscretizationLocal.STEPS) ;
nparam = length(INPUT_PARAMETERS{1})-1 ;
MATRIX_POINTS_SPACE_PARAMETER = zeros(nstepsLOC*length(CASES),nparam) ;
TESTING_PARAM = [] ; 
INDEXES_TRAJECTORIES = {} ; 
iacum_proj = 1 ;
for iproj = 1:length(CASES)    
    disp(['iproj = ',num2str(CASES(iproj))])  
    
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    if CASES(iproj) == 1
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'MATPRO','MESH','OPERFE','DISP_CONDITIONS',...
            'OTHER_output') ;
        DATA = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ;
        DOFl = DISP_CONDITIONS.DOFl ;
        DOFr = DISP_CONDITIONS.DOFr ;
    end    
    TIME_STEPS{iproj} = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STEPS ;    
    
    % INPUT PARAMETER
    % **************
    MU_1 = INPUT_PARAMETERS{iproj}(1)*TIME_STEPS{iproj}   ;
    MU_2 = INPUT_PARAMETERS{iproj}(2)*TIME_STEPS{iproj}  ;
    
     
    if ~isempty(DATAoffline.TESTING_TRAJ_IS_TRAINING_TRAJ)
       if  DATAoffline.TESTING_TRAJ_IS_TRAINING_TRAJ == CASES(iproj)
           TESTING_PARAM = [MU_1',MU_2'] ; 
       end
    end 
    
    iacum_proj_new = iacum_proj + length(MU_1)-1;
    INDEXES_TRAJECTORIES{iproj} = iacum_proj:iacum_proj_new ; 
    MATRIX_POINTS_SPACE_PARAMETER(iacum_proj:iacum_proj_new,:) =[MU_1',MU_2'] ;
    iacum_proj =   iacum_proj_new + 1;
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)        
        %DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;        
        % Or the whole matrix ....
        disp(['Rank of the matrix = ',num2str(length(SNAP_cluster.DISP.S))]) ;
        
%         if DATAoffline.NormalizedTrainingSnapshots == 1
%             SNAP_cluster.DISP.S = SNAP_cluster.DISP.S/SNAP_cluster.DISP.S(1) ; 
%         end
        
        DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        DISP_LOC{iloc} =  DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;   
        
        
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);    
    %      if    DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
    %          SNAPdisp_ALL{iproj} =SNAPdisp{iproj} ;
    %      end    
   % if DATAoffline.PLOT_MODES_WITH_DOFR == 0 & DATAoffline.CLUSTER.WITH_ALL_DOFS==0
    %    SNAPdisp{iproj} = SNAPdisp{iproj}(DOFl,:) ;        
   % end    
end
% Now we apply the partioned SVD
% ----------------------------
TOL_BLOCK = DATAoffline.TOL_SVD_MATRIX_ALL*ones(length(SNAPdisp),1)' ;

%   if DATAoffline.PLOT_MODES_WITH_DOFR == 0 && DATAoffline.CLUSTER.WITH_ALL_DOFS==0
%       DOFsINCLUDE = DOFl ; 
%   else
%       DOFsINCLUDE
%   end

DATAsvd=[];
DATAsvd.HIDE_OUTPUT = 1 ;

% This is the basis matrix for the case of just one cluster ---DONE WITH
% ALL THE DOFS !!!!!

[BasisU,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPdisp,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of displacement modes =',num2str(size(BasisU,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])