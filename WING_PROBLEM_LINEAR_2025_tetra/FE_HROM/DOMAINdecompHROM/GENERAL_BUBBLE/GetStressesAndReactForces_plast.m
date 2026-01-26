function [SNAPreactFORCESproj,SNAPstressSTWOproj,SNAPstressPonePROJ,NAME_BASE]  =...
    GetStressesAndReactForces_plast(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
    INFO_RVE,BasisDEFdisp,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB)
% Copy of GetStressesAndReactForces_1dom.m
% Adapted to cope with "bubble" modes (complementary modes)
if nargin == 0
    load('tmp.mat')
end

% CHECK_ERROR_NONBUBBLE = 0;
%
% if CHECK_ERROR_NONBUBBLE == 1
%    disp('---------------------------------------------------------------------------')
%    disp(['Disp. matrix formed by non-bubble functions (orthogonal complement). Checking stress errors'])
%    disp('---------------------------------------------------------------------------')
%     [BasisU,~,~] = SVDT([PhiRB,BasisDEFdisp.BASIC,BasisDEFdisp.BasisDEFdispcomp_orth]) ;
%
% else
%      disp('---------------------------------------------------------------------------')
%    disp(['Disp. matrix formed by RB + BASIC. DEF + Bubble modes.  '])
%    disp('---------------------------------------------------------------------------')
DATAoffline = DefaultField(DATAoffline,'CHECK_ACCURACY_PROJECTED_STRESSES',1) ;

if DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES ==1
    [BasisU,~,~] = SVDT([PhiRB,BasisDEFdisp.UpsilonDEF]) ;
    
else
    disp(['Stresses and reactive forces arising from exact displacements (error must be zero)'])
    BasisU = 1; 
end
%end

% ------------------------------------------------------------------------------------------------------
% Now we have to compute the PK2 stresses produced by these strains, and check
% that all blocks are approximated to an accuracy theshold
% DATAoffline.\errorSTRESS
CASES = 1:size(DATAcommon.INPUTS_PARAMETERS,2) ;

STRESS_PK2_error = zeros(length(CASES),1) ;

SNAPstressSTWOproj = cell(1,length(CASES)) ;
SNAPstressPonePROJ = cell(1,length(CASES)) ;
SNAPreactFORCESproj = cell(1,length(CASES)) ;

VAR = [] ; VARint_n = [] ;
DATA = DefaultField(DATA,'ListFieldInternalVariables',[] ) ;



 

DATAoffline = DefaultField(DATAoffline,'AdditionalTests',[] );

if ~isempty(DATAoffline.AdditionalTests) ;
    CASES = 1:(size(DATAcommon.INPUTS_PARAMETERS,2) + length(DATAoffline.AdditionalTests)) ;
end


for iproj = 1:length(CASES)
    % NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    
    if iproj <= size(DATAcommon.INPUTS_PARAMETERS,2)
        NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
        ISADD = 0 ;
    else
        iloc = iproj-size(DATAcommon.INPUTS_PARAMETERS,2) ;
        NAME_FOLDER = DATAoffline.AdditionalTests{iloc} ;
        ISADD = 1;
    end
    
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    
    
    
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    SNAPstressSTWOproj_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SNAPstressSTWO_LOC  = cell(1,length(NAME_SNAP_loc)) ;
    
    
    
    SNAPstressPonePROJ_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SELECTED_ENTRIES_GAUSS = INFO_RVE.GaussINDEX_stress ;
    SELECTED_ENTRIES_DOFS = INFO_RVE.DOFS_globNUM ;
    
    if ISADD ==0 && ~isempty(Fbody)
        Fbody_loc = Fbody{iproj} ;
        Ftrac_loc = Ftrac{iproj} ;
    else
        % Additional training test
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'OTHER_output') ;
        Fbody_loc =OTHER_output.Fbody ;
        Ftrac_loc =OTHER_output.Ftrac ;
        if ~any(Fbody_loc.U)
            Fbody_loc.U = [] ;
            Fbody_loc.a = [] ;
        end
        
    end
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        
        % 2ND PK STRESSES (Cauchy stresses for small strains )
        SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U(SELECTED_ENTRIES_GAUSS,:)',SNAP_cluster.PK2STRESS.S)' ;
        SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        % Projected displacements
        % -----------------------
      %  if ~isempty(BasisU)
        coeff = BasisU'*SNAP_cluster.DISP.U(SELECTED_ENTRIES_DOFS,:) ;
        coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
        coeff = coeff*SNAP_cluster.DISP.V' ;
        % ------------------------------------------------
        d = BasisU*coeff; % Snapshot displacements (local)
      %  else
            
       % end
        
        % Computation of stresses and reactive forces (residual) arising
        % from the projected stresses. Valid for both linear and  nonlinear
        % regimes
        [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR,SNAPreactFORCESproj{iproj}] = MultiSnapStressFromDispLOC...
            (VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
            SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj,INFO_RVE,Fbody_loc,Ftrac_loc) ;
        
    end
    SNAPstressSTWOproj_LOC = cell2mat(SNAPstressSTWOproj_LOC) ;   % Approximate
    SNAPstressSTWO_LOC = cell2mat(SNAPstressSTWO_LOC) ;   % exact (from FE training tests)
    
    % Recall that we may not want to study the whole range of snapshots ()
%     nsnapshots = ceil(size(SNAPstressSTWO_LOC,2)*min(DATAoffline.proportionSNAPSHOTS(iproj),1)) ;
%     SNAPstressSTWO_LOC = SNAPstressSTWO_LOC(:,1:nsnapshots) ;
%     SNAPstressSTWOproj_LOC = SNAPstressSTWOproj_LOC(:,1:nsnapshots) ;
%     
    % Check if it meets the error criterion
    STRESS_PK2_error(iproj) = norm(SNAPstressSTWO_LOC-SNAPstressSTWOproj_LOC,'fro')/norm(SNAPstressSTWO_LOC,'fro')  ;
    disp(['Project = ',num2str(iproj),'; ERROR stress PK2= ',num2str(STRESS_PK2_error(iproj))]);
    
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
  %  [UU,SS,VV] =  RSVDT(SNAPreactFORCESproj{iproj}(:,1:nsnapshots)) ;
  %s  SNAPreactFORCESproj{iproj} = bsxfun(@times,UU',SS/SS(1))' ;
    
 %   if DATAoffline.BASIS_STRESSES_FROM_PROJECTED_DISPLACEMENTS ==1
        % THESE ARE THE STRESSES AND THE INTERACTION FORCES ARISING FROM
        % THE PROJECTED DISPLACEMENTS
       % [UU,SS,VV] =  RSVDT(SNAPstressSTWOproj_LOC) ;
       % SNAPstressSTWOproj{iproj} = bsxfun(@times,UU',SS/SS(1))' ;  % This is equivalent to NORMALIZATION
        if DATAoffline.BASIS_STRESSES_FROM_PROJECTED_DISPLACEMENTS ==1
         SNAPstressSTWOproj{iproj} = SNAPstressSTWOproj_LOC ;  
        else
             SNAPstressSTWOproj{iproj} = SNAPstressSTWO_LOC ;  
        end
        
  %  else
        % tHESE ARE THE "EXACT" VALUES OF STRESSES (TRAINING TESTS)
   %     [UU,SS,VV] =  RSVDT(SNAPstressSTWO_LOC) ;
    %    SNAPstressSTWOproj{iproj} = bsxfun(@times,UU',SS/SS(1))' ;  % This is equivalent to NORMALIZATION
        
   % end
    
    
end

if  DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES ==1
   warning(['RUN AGAIN  SETTING  DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES = 0']) ;  
   error('')
 
 end

%
% if CHECK_ERROR_NONBUBBLE == 1
%    error('Set  CHECK_ERROR_NONBUBBLE = 0 to proceed' )
% end


% % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK2 STRESSES (FOR RECONSTRUCTION PURPOSES)
% % -------------------------------------------------------------------------------------------------
% if length(DATAoffline.errorPK2stress_basis) == 1
%     TOL_BLOCK = DATAoffline.errorPK2stress_basis*ones(length(SNAPstressSTWOproj),1)' ;
% else
%     TOL_BLOCK = DATAoffline.errorPK2stress_basis ;
% end
%
% USE_BLOCKED_VERSIONRandomized =0 ;
%
% if USE_BLOCKED_VERSIONRandomized == 1
%     DATAsvd=[];
%     [BasisStwo,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressSTWOproj,TOL_BLOCK,DATAsvd) ;
%     disp('***********************************************************')
%     disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
%     disp('***********************************************************')
%
%
%     % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK1 STRESSES (FOR applying the ECM )
%     % ----------------------------------------------------------------------------------------------
%     DATAsvd=[];
%     [BasisPone,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressPonePROJ,TOL_BLOCK,DATAsvd) ;
%     % disp('***********************************************************')
%     % disp(['Number of PK1 stress modes =',num2str(size(BasisPone,2))])
%     % disp('***********************************************************')
%     % The basis for the PK1 stresses is weighted with the singular values
%
%
% else
%
%     % NORMALIZED STRESSES...
%     SNAPstressSTWOproj = cell2mat(SNAPstressSTWOproj) ;
%     nSNAPstressSTWOproj = sqrt(sum(SNAPstressSTWOproj.^2,1));
%     SNAPstressSTWOproj = bsxfun(@times,SNAPstressSTWOproj',1./nSNAPstressSTWOproj')' ;
%
%     TOL_BLOCK = DATAoffline.errorPK2stress_basis  ;
%     DATAsvd.RELATIVE_SVD=1;
%     [BasisStwo,S,V] = SVDT((SNAPstressSTWOproj),TOL_BLOCK,DATAsvd) ;
%     disp('***********************************************************')
%     disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
%     disp('***********************************************************')
%
%     BasisPone = BasisStwo ;
%
%     % [BasisPone,S,V] = SVDT(cell2mat(SNAPstressPonePROJ),TOL_BLOCK,DATAsvd) ;
%     %     disp('***********************************************************')
%     %     disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
%     %     disp('***********************************************************')
%     % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK1 STRESSES (FOR applying the ECM )
%     % ----------------------------------------------------------------------------------------------
%
%     % [BasisPone,S,V] = SVDT(SNAPstressPonePROJ,TOL_BLOCK,DATAsvd) ;
%     % disp('***********************************************************')
%
% end
% DATAoffline = DefaultField(DATAoffline,'IncludeSingularValueBasisStressesPone',1) ; % =
% if DATAoffline.IncludeSingularValueBasisStressesPone == 1
%     BasisPone = bsxfun(@times,BasisPone',S)' ;
% end
%
%
%
%
%
% % NORMALIZATION
% for iproj = 1:length(SNAPreactFORCESproj)
%     SNAPreactFORCESproj{iproj} = SNAPreactFORCESproj{iproj}/norm(SNAPreactFORCESproj{iproj},'fro') ;
% end
% SNAPreactFORCESproj = cell2mat(SNAPreactFORCESproj) ;
%
%
% % REACTIVE FOCES, BASIS MATRICES (not strictly necessary)
% % ------------------------------
% % DATAoffline = DefaultField(DATAoffline,'errorREACTIVE_FORCES',0) ;
% % TOL_BLOCK = DATAoffline.errorREACTIVE_FORCES*ones(length(SNAPreactFORCESproj),1)' ;
% % THIS STEP IS NOT STRICLTLY NECESSARY !!!!
%
% DATAoffline= DefaultField(DATAoffline,'PlotRawReactions',0); % = 1 ;
%
%
%
% if DATAoffline.PlotRawReactions == 1
%
%     % [BasisREACTF,S,V] = SVDT(SNAPreactFORCESproj) ;
%
%     BasisREACTF = SNAPreactFORCESproj(:,2:2:end) ;
%
%     NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
%     if  ~exist(NAME_MODES_FOLDER)
%         mkdir(NAME_MODES_FOLDER)
%     end
%     NAME_MODES_DISP = [NAME_MODES_FOLDER,NAME_BASE,'REACTmodesRAW'] ;
%
%     NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
%     NameFile_res = [NAME_MODES_DISP,'.res'] ;
%     DATALOC = [] ;
%     GidPostProcessModesDOML(MESHdom.COOR,MESHdom.CN,MESHdom.TypeElement,BasisREACTF,DATA.MESH.posgp,NameFileMesh,NameFile_res,MESHdom.MaterialType,DATALOC) ;
%
% end