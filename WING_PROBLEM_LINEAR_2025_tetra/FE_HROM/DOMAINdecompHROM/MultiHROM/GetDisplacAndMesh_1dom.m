function [INFO_RVE,MESHdom,SNAPdisp,BasisU,DATA,OTHER_output,OPERFE,MATPRO,Fbody,Ftrac] ...
    = GetDisplacAndMesh_1dom(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline)
% This function returns the finite element information (mesh and displacements) stored in data
% structure INFO_SNAPSHOTS of a given parametric finite element
% information.
% VERSION 8-Feb-2023:  It returns DISPLACEMENTS MODES and MESH information
% of one single subdomain
% -----------------------------------
if nargin == 0
    load('tmp.mat')
end

Fbody = [] ; Ftrac = [] ; 

CASES = 1:size(DATAcommon.INPUTS_PARAMETERS,2) ;

DATAoffline = DefaultField(DATAoffline,'AdditionalTests',[] );
DATAoffline = DefaultField(DATAoffline,'ReturnFullDisplacementMatrix',0 );

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
        DATAoffline = DefaultField(DATAoffline,'LABEL_DOMAIN',1) ; 
        IDOM = DATAoffline.LABEL_DOMAIN; % Domain we wish to study
        DATAoffline = DefaultField(DATAoffline,'LABELS_FACES',1:length(MESH.NODES_FACES)) ;  % 17-Nov-2023 (not verified properly...)
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
        
    %    INCLUDE_V = 0;
        if DATAoffline.ReturnFullDisplacementMatrix == 1
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


% THIS STEP IS NOT STRICTLY NECESSARY
DATAoffline = DefaultField(DATAoffline,'USE_BASIS_MATRIX_DISPLACEMENTS_AS_SNAPSHOTS',0) ; % = 1;

if DATAoffline.USE_BASIS_MATRIX_DISPLACEMENTS_AS_SNAPSHOTS == 1

DATAsvd=[];
[BasisU,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPdisp,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of displacement modes =',num2str(size(BasisU,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])
else
    BasisU = 1 ; U = [] ; S = [] ; 
end



if DATAoffline.USE_BASIS_MATRIX_DISPLACEMENTS_AS_SNAPSHOTS == 0
    SNAPdisp = cell2mat(SNAPdisp) ;
elseif DATAoffline.USE_BASIS_MATRIX_DISPLACEMENTS_AS_SNAPSHOTS == 1
    SNAPdisp = BasisU ;

end



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
