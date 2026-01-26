function [INFO_RVE,MESHdom,SNAPdisp,BasisU,DATA,OTHER_output,OPERFE,MATPRO,Fbody,Ftrac] ...
    = GetDisplacAndMesh_1domBUB(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline)
% This function returns the finite element information (mesh and displacements) stored in data
% structure INFO_SNAPSHOTS of a given parametric finite element
% information.
% VERSION 8-Feb-2023:  It returns DISPLACEMENTS MODES and MESH information
% of one single subdomain
% Bubble version, June-05-2023
% -----------------------------------
if nargin == 0
    load('tmp1.mat')
end


CASES = 1:size(DATAcommon.INPUTS_PARAMETERS,2) ;

DATAoffline = DefaultField(DATAoffline,'AdditionalTests',[] );

if ~isempty(DATAoffline.AdditionalTests) ;
    CASES = 1:(size(DATAcommon.INPUTS_PARAMETERS,2) + length(DATAoffline.AdditionalTests)) ;
end



% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
for iproj = 1:length(CASES)
    
    if iproj <= size(DATAcommon.INPUTS_PARAMETERS,2)
        NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    else
        iloc = iproj-size(DATAcommon.INPUTS_PARAMETERS,2) ;
        NAME_FOLDER = DATAoffline.AdditionalTests{iloc} ;
    end
    
    
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    if iproj == 1
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'MATPRO','MESH','OPERFE','OTHER_output','Fbody','Ftrac') ;
        %
        
        DATA = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ;
        
        % --------------------------------------------------------
        IDOM = DATAoffline.LABEL_DOMAIN; % Domain we wish to study
        IndexBoundaries = DATAoffline.LABELS_FACES ;% Indexes boundaries corresponding to this domain
        [INFO_RVE,MESHdom,OPERFE,MATPRO] = IndexesREFdomain(IDOM,IndexBoundaries,MESH,DATA,OPERFE,MATPRO) ; % Indexes used for extracting information
    end
    % load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'OTHER_output')
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SELECTED_DOFS = INFO_RVE.DOFS_globNUM ;
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
        DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U(SELECTED_DOFS,:)',SNAP_cluster.DISP.S)' ;
        
        INCLUDE_V = 0;
        if INCLUDE_V == 1
            DISP_LOC{iloc} = DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
        end
        %
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    
end


%


% Now we apply the partioned SVD
% ----------------------------
if length(DATAoffline.errorDISP) == 1
    TOL_BLOCK = DATAoffline.errorDISP*ones(length(SNAPdisp),1)' ;
else
    TOL_BLOCK = DATAoffline.errorDISP ;
end



% NORMALIZATION
for iproj = 1:length(SNAPdisp)
    SNAPdisp{iproj} = SNAPdisp{iproj}/norm(SNAPdisp{iproj},'fro') ; %
end

%%% NOW WE HAVE ELASTIC SNAPSHOTS AND BUBBLE SNAPSHOTS
% ----------------------------------------------------
% From 1 to end-length(DATAoffline.AdditionalTests): elastic "basic"
% snapshots
% Then ---> Bubble snapshots 

% BASIC SNAPSHOTS 
%--------------------------------------------------
INDbasicSNAP = 1:(length(CASES)-length(DATAoffline.AdditionalTests)); 
TOL_BLOCK_basic = TOL_BLOCK(INDbasicSNAP);

DATAsvd=[];
[BasisU_basic,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPdisp(INDbasicSNAP),TOL_BLOCK_basic,DATAsvd) ;
disp('***********************************************************')
disp(['Number of displacement modes (basic) =',num2str(size(BasisU_basic,2)) ])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])

% BUBBLE MODES 
% ------------
INDBUBBLE =  (length(CASES)-length(DATAoffline.AdditionalTests)+1):length(CASES); 
TOL_BLOCK_BUBLE= TOL_BLOCK(INDBUBBLE);

[BasisU_BUBBLE,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPdisp(INDBUBBLE),TOL_BLOCK_BUBLE,DATAsvd) ;
disp('***********************************************************')
disp(['Number of displacement modes (BUBBLE) =',num2str(size(BasisU_BUBBLE,2))])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])


%%%%%%% 
% INTERSECTION BETWEEN BASISU_BUBBLE AND BasisU_basic
[uu,ss,vv] = SVDT(BasisU_basic'*BasisU_BUBBLE) ; 

disp(['Cosine Angles subtended (in degrees) by basis functions for buble and basic snapshots '])
ss

disp('Using basis modes for basic snapshots as basis matrix for all snapshots')
BasisU = BasisU_basic ; 

% Snapshots 
SNAPdisp = cell2mat(SNAPdisp(INDbasicSNAP)) ;  % Only basic modes are included here 



% ***************************
% PLOTTING DISPLACEMENT MODES
% ***************************

NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
if  ~exist(NAME_MODES_FOLDER)
    mkdir(NAME_MODES_FOLDER)
end

PLOT_RAW_MODES = 0 ;

if PLOT_RAW_MODES == 1
    NAME_MODES_DISP = [NAME_MODES_FOLDER,NAME_BASE,'DISPmodesRAW'] ;
    
    NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
    NameFile_res = [NAME_MODES_DISP,'.res'] ;
    DATALOC = [] ;
  %  GidPostProcessModesDOML(MESHdom.COOR,MESHdom.CN,MESHdom.TypeElement,SNAPdisp(:,2:2:end),DATA.MESH.posgp,NameFileMesh,NameFile_res,MESHdom.MaterialType,DATALOC) ;
    GidPostProcessModesDOML(MESHdom.COOR,MESHdom.CN,MESHdom.TypeElement,SNAPdisp ,DATA.MESH.posgp,NameFileMesh,NameFile_res,MESHdom.MaterialType,DATALOC) ;

end
