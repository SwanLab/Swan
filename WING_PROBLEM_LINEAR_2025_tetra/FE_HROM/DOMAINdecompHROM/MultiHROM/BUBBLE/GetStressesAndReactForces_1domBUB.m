function [SNAPreactFORCESproj,BasisStwo,BasisPone,NAME_BASE]  =...
    GetStressesAndReactForces_1domBUB(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
    INFO_RVE,MESHdom,BasisU,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output)
% THIS FUNCTIO READS THE STRESS INFORMATION CORRESPONDING TO THE TRAINING
% PROJECTS, AND RETRIEVE THE STRESSES CORRESPONDING TO THE SELECTED
% SUBDOMAIN
% IT ALSO PROJECTS THE DISPLACEMENTS HISTORY OF EACH PROJECT ONTO THE BASIS
% MATRIX BasisU, and compute the difference between the FE stresses and the
% one determined from the approximated displacement history
% Lastly, it also calculates a  basis matrix for the column space of the
% reactive snapshots (in turn, such snapshots are determined as REACTIVE_forces = InternalForces-ExternalForces)
% Joaquín A. Hernández, 9-Feb-2023
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx

% Bubble version: /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/09_POSTPAPER_POROUS/02_bubbleMODES.mlx

% ------------------------------------------------------------------------------------------------------
% Now we have to compute the PK2 stresses produced by these strains, and check
% that all blocks are approximated to an accuracy theshold
% DATAoffline.\errorSTRESS

if nargin == 0 
    load('tmp1.mat')
end 

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
    
    if ISADD ==0
        Fbody_loc = Fbody{iproj} ;
        Ftrac_loc = Ftrac{iproj} ;
    else
        % Additional training test
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'OTHER_output') ;
        Fbody_loc =OTHER_output.Fbody ;
        Ftrac_loc =OTHER_output.Ftrac ;
    end
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        
        % 2ND PK STRESSES (Cauchy stresses for small strains )
        SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U(SELECTED_ENTRIES_GAUSS,:)',SNAP_cluster.PK2STRESS.S)' ;
        SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        % Projected displacements
        % -----------------------
        coeff = BasisU'*SNAP_cluster.DISP.U(SELECTED_ENTRIES_DOFS,:) ;
        coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
        coeff = coeff*SNAP_cluster.DISP.V' ;
        % ------------------------------------------------
        d = BasisU*coeff; % Snapshot displacements (local)
        
        %         d_noproyec = SNAP_cluster.DISP.U(SELECTED_ENTRIES_DOFS,:)  ;
        %         coeff =bsxfun(@times,d_noproyec',SNAP_cluster.DISP.S)' ;
        %         d_noproyec = coeff*SNAP_cluster.DISP.V' ;
        %         nddd = norm(d_noproyec-d,'fro') ;
        
        %         dR = bsxfun(@times,SNAP_cluster.DISP.U(DOFr,:)',SNAP_cluster.DISP.S)' ;
        %         dR = dR*SNAP_cluster.DISP.V' ;
        %         ndof = size(dL,1)+size(dR,1) ;
        %         d = zeros(ndof,size(dL,2)) ;
        %         d(DOFl,:)  = dL ;
        %         d(DOFr,:)  = dR ;
        %
        % 2. Deformation gradient at all Gauss points
        % SMALL STRAINS
        
        
        
        [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR,SNAPreactFORCESproj{iproj}] = MultiSnapStressFromDispLOC...
            (VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
            SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj,INFO_RVE,Fbody_loc,Ftrac_loc) ;
        
    end
    SNAPstressSTWOproj_LOC = cell2mat(SNAPstressSTWOproj_LOC) ;   % Approximate
    SNAPstressSTWO_LOC = cell2mat(SNAPstressSTWO_LOC) ;   % exact
    
    % Check if it meets the error criterion
    STRESS_PK2_error(iproj) = norm(SNAPstressSTWO_LOC-SNAPstressSTWOproj_LOC,'fro')/norm(SNAPstressSTWO_LOC,'fro')  ;
    disp(['Project = ',num2str(iproj),'; ERROR stress PK2= ',num2str(STRESS_PK2_error(iproj))]);
    %     if STRESS_PK2_error(iproj) > DATAoffline.errorSTRESS
    %         %  dbstop('129')
    %         error('STRESS ERROR CRITERION NOT MET: TRY WITH A LOWER ERROR TOLERANCE FOR DISPLACEMENTS ')
    %     end
    
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
    % Change 6-jun-2023 
    % SNAPstressSTWOproj_LOC --> Approximated 
    % SNAPstressSTWO_LOC --> Exact 
  %  [UU,SS,VV] =  RSVDT(SNAPstressSTWOproj_LOC) ;
 
        [UU,SS,VV] =  RSVDT(SNAPstressSTWO_LOC) ;
 
    SNAPstressSTWOproj{iproj} = bsxfun(@times,UU',SS)' ;
    % PK1 stresses
%     % ----------------
%     SNAPstressPonePROJ_LOC = cell2mat(SNAPstressPonePROJ_LOC) ;   % exact
%     [UU,SS,VV] =  RSVDT(SNAPstressPonePROJ_LOC) ;
    SNAPstressPonePROJ{iproj} = SNAPstressSTWOproj{iproj} ; % Makes no sense to maintain this distinction in the small strain regime 
    
    
end

% BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK2 STRESSES (FOR RECONSTRUCTION PURPOSES)
% -------------------------------------------------------------------------------------------------
if length(DATAoffline.errorPK2stress_basis) == 1
    TOL_BLOCK = DATAoffline.errorPK2stress_basis*ones(length(SNAPstressSTWOproj),1)' ;
else
    TOL_BLOCK = DATAoffline.errorPK2stress_basis ;
end

USE_BLOCKED_VERSIONRandomized =0 ;

if USE_BLOCKED_VERSIONRandomized == 1
    DATAsvd=[];
    [BasisStwo,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressSTWOproj,TOL_BLOCK,DATAsvd) ;
    disp('***********************************************************')
    disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
    disp('***********************************************************')
    
    
    % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK1 STRESSES (FOR applying the ECM )
    % ----------------------------------------------------------------------------------------------
    DATAsvd=[];
    [BasisPone,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressPonePROJ,TOL_BLOCK,DATAsvd) ;
    % disp('***********************************************************')
    % disp(['Number of PK1 stress modes =',num2str(size(BasisPone,2))])
    % disp('***********************************************************')
    % The basis for the PK1 stresses is weighted with the singular values
    
    
else
    
    % NORMALIZED STRESSES...
    SNAPstressSTWOproj = cell2mat(SNAPstressSTWOproj) ;
    nSNAPstressSTWOproj = sqrt(sum(SNAPstressSTWOproj.^2,1));
    SNAPstressSTWOproj = bsxfun(@times,SNAPstressSTWOproj',1./nSNAPstressSTWOproj')' ;
    
    TOL_BLOCK = DATAoffline.errorPK2stress_basis  ;
    DATAsvd.RELATIVE_SVD=1;
    [BasisStwo,S,V] = SVDT((SNAPstressSTWOproj),TOL_BLOCK,DATAsvd) ;
    disp('***********************************************************')
    disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
    disp('***********************************************************')
    
    BasisPone = BasisStwo ;
    
    % [BasisPone,S,V] = SVDT(cell2mat(SNAPstressPonePROJ),TOL_BLOCK,DATAsvd) ;
    %     disp('***********************************************************')
    %     disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
    %     disp('***********************************************************')
    % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK1 STRESSES (FOR applying the ECM )
    % ----------------------------------------------------------------------------------------------
    
    % [BasisPone,S,V] = SVDT(SNAPstressPonePROJ,TOL_BLOCK,DATAsvd) ;
    % disp('***********************************************************')
    
end
DATAoffline = DefaultField(DATAoffline,'IncludeSingularValueBasisStressesPone',1) ; % =
if DATAoffline.IncludeSingularValueBasisStressesPone == 1
    BasisPone = bsxfun(@times,BasisPone',S)' ;
end

INDbasicSNAP = 1:(length(CASES)-length(DATAoffline.AdditionalTests)); %Indices basic snapshots
INDBUBBLE =  (length(CASES)-length(DATAoffline.AdditionalTests)+1):length(CASES); 




% NORMALIZATION ---REACTIVE MODES 
for iproj = 1:length(SNAPreactFORCESproj)
    SNAPreactFORCESproj{iproj} = SNAPreactFORCESproj{iproj}/norm(SNAPreactFORCESproj{iproj},'fro') ;
end

disp(['Only basic reactive modes are considered'])
SNAPreactFORCESproj = cell2mat(SNAPreactFORCESproj(INDbasicSNAP)) ;  % Only reactive basic modes are considered


% REACTIVE FOCES, BASIS MATRICES (not strictly necessary)
% ------------------------------
% DATAoffline = DefaultField(DATAoffline,'errorREACTIVE_FORCES',0) ;
% TOL_BLOCK = DATAoffline.errorREACTIVE_FORCES*ones(length(SNAPreactFORCESproj),1)' ;
% THIS STEP IS NOT STRICLTLY NECESSARY !!!!

PLOTmodes = 0 ;

if PLOTmodes == 1
    
    [BasisREACTF,S,V] = SVDT(SNAPreactFORCESproj) ;
    
    % BasisREACTF = SNAPreactFORCESproj(:,2:2:end) ;
    
    NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
    if  ~exist(NAME_MODES_FOLDER)
        mkdir(NAME_MODES_FOLDER)
    end
    NAME_MODES_DISP = [NAME_MODES_FOLDER,NAME_BASE,'REACTmodesRAW'] ;
    
    NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
    NameFile_res = [NAME_MODES_DISP,'.res'] ;
    DATALOC = [] ;
    GidPostProcessModesDOML(MESHdom.COOR,MESHdom.CN,MESHdom.TypeElement,BasisREACTF,DATA.MESH.posgp,NameFileMesh,NameFile_res,MESHdom.MaterialType,DATALOC) ;
    
end