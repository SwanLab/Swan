function [DATAIN] = PrintingGid_ECMpoints(DATAIN,DATA_REFMESH,HYPERREDUCED_VARIABLES,HISTORY_ITERATIONS) ;

if nargin == 0
    load('tmp1.mat')
    
end

% ---------------------------------------

CNref = DATA_REFMESH.CN  ;
COOR = DATA_REFMESH.COOR  ;

DATA= [] ;


% Name of the mesh fileRpo
[DIRinp,NAMEFILEinput ]=fileparts(DATAIN.NAME_WS_MODES) ;
DATAIN = DefaultField(DATAIN,'LABEL_NAME_PROJECT','') ;
NAMEFILEinput=  ['ECM_ELEMENTS',DATAIN.LABEL_NAME_PROJECT] ;

ROOTFOLDER = [DIRinp,filesep,'GIDPOST',filesep] ;
disp(['Writing ECM mesh  in folder: ',ROOTFOLDER]);
if ~exist(ROOTFOLDER)
    mkdir(ROOTFOLDER) ;
end
NameFile_msh = [ROOTFOLDER,filesep,NAMEFILEinput,'' ,'.msh'] ;
NameFile_res = [ROOTFOLDER,filesep,NAMEFILEinput,'' ,'.res'] ;


COORprint = COOR ;
CNprint = {CNref} ;
NAME_INPUT_DATA =NAMEFILEinput ;
MaterialTypegloPRINT = {DATA_REFMESH.MaterialType} ;
TypeElementPRINT = {DATA_REFMESH.TypeElement} ;
NAMEMESH = {'FULL'} ;


IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COORprint,CNprint,NAME_INPUT_DATA,MaterialTypegloPRINT,...
    TypeElementPRINT,NAMEMESH);
%    WEIGHTS = HYPERREDUCED_VARIABLES.WdomRED ;


DATA.NAMEMESH = NAMEMESH ;
DATA.WEIGHTS = HISTORY_ITERATIONS.WEIGHTS;
DATA.ELEMENTS_WITH_WEIGHTS = HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS;

posgp = {[],[]};
DATAINloc = [] ;

[NameFileRES ]= GidPostProcess_WEIGHTSHISTORY(COORprint,CNprint,TypeElementPRINT, ...
    posgp,NameFile_res,MaterialTypegloPRINT,IND_ELEM_MESHES,...
    DATA_REFMESH,DATA,DATAINloc);




disp('-----------------------------')
disp(['To see reduced set of integration points, open']);
disp(NameFile_msh);
disp(['--------------------- '])

%
%     GidPostProcess_onlyECMpoints(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdef,DATA_REFMESH.posgp,...
%         NAME_MODES,DATA,LEGENDG);
DATAIN = DefaultField(DATAIN,'MSGPRINT',{});
DATAIN.MSGPRINT{end+1} = '---------------------------------------------------' ;
DATAIN.MSGPRINT{end+1} = ['To see reduced set of integration points, open'] ;
DATAIN.MSGPRINT{end+1} =NameFile_msh ;


