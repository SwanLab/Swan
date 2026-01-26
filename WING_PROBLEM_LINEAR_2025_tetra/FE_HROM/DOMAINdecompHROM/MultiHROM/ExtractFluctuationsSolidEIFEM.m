function dFLUCT = ExtractFluctuationsSolidEIFEM(INPUT_PROBLEMS,DATAoffline,ILOADSTATES,...
    MESH,MATPRO,OPERFE,DISP_CONDITIONS,OTHER_output) ;

if nargin == 0
    load('tmp.mat')
end
%
% READING INPUT DATA
DATAcommon = feval(INPUT_PROBLEMS) ;
DATAcommon = DefaultField(DATAcommon,'NameFileMeshDATA',[] ) ;
DATAoffline = DefaultField(DATAoffline,'NAMEMESH', DATAcommon.NameFileMeshDATA);
DATAcommon = DefaultField(DATAcommon,'MESH_QUADRATIC_PARENT_DOMAIN', []);
% JAHO, 30-Nov-2023
DATAoffline = DefaultField(DATAoffline,'MESH_QUADRATIC_PARENT_DOMAIN',DATAcommon.MESH_QUADRATIC_PARENT_DOMAIN) ;
DATAcommon.NameFileMeshDATA = DATAoffline.NAMEMESH ;
NAME_BASE =[DATAcommon.NameParamStudyLOC,'_param_'];
NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];
% Let us construct first the matrix the snapshots of each project
SNAPdisp =[] ;
iproj = ILOADSTATES ;
NAME_FOLDER = [NAMEsnap_base,num2str(iproj)] ;
NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
load(NAME_INFO,'INFO_SNAPSHOTS')
%load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'MATPRO','MESH','OPERFE','OTHER_output','Fbody','Ftrac') ;
%DATA = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ;
% --------------------------------------------------------
%DATAoffline = DefaultField(DATAoffline,'LABEL_DOMAIN',1) ;
%IDOM = DATAoffline.LABEL_DOMAIN; % Domain we wish to study
%DATAoffline = DefaultField(DATAoffline,'LABELS_FACES',1:length(MESH.NODES_FACES)) ;  % 17-Nov-2023 (not verified properly...)
% load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'OTHER_output')

STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
 for iloc = 1:length(NAME_SNAP_loc)
    Nameloc = NAME_SNAP_loc{iloc} ;
    load(Nameloc,'SNAP_cluster') ;
    % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
    DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        DISP_LOC{iloc} = DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
 end

 % DISPLACEMENT VECTOR 
d =   DISP_LOC{1}(:,2);

% FLUCTUATIONS 
bDOFS = MESH.INFO_PERIODIC_CONDITIONS.boundaryDOFS ;
dFLUCT = d(bDOFS)-MESH.INFO_PERIODIC_CONDITIONS.ImposedDisplacement  ;

ndim = size(MESH.COOR,2) ; 

BoundaryNodes = large2small(bDOFS,ndim) ; 
COORbnd = MESH.COOR(BoundaryNodes,:) ;

CNb = cell(length(MESH.Indexes_faces_bnd_element),1) ;
for iface  = 1:length(MESH.Indexes_faces_bnd_element)
    CNb{iface} = MESH.CNb(MESH.Indexes_faces_bnd_element{iface},:) ;
end
CNb = cell2mat(CNb(:)) ;
CNbREN  =  RenumberConnectivities( CNb,1:length(BoundaryNodes) );


NameLoc =     ['IMPOSED_DISPT_AND_FLUCT_',num2str(ILOADSTATES)] ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = []; 

 GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[MESH.INFO_PERIODIC_CONDITIONS.ImposedDisplacement,d(bDOFS),dFLUCT],...
[],NameFileMesh,NameFile_res,[],DATALOC);


