function PlotModesStress_SUBDOMAIN(DATA,BasisSTRESS,MESH)

if nargin == 0
    load('tmp1.mat')
end

NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
if  ~exist(NAME_MODES_FOLDER)
    mkdir(NAME_MODES_FOLDER)
end

DATA = DefaultField(DATA,'Basic','LEGEND_MODES_stress') ; 

NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,'STRESSmodes',DATA.LEGEND_MODES_stress] ;

NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [] ;

DATA.MaterialType = MESH.MaterialType ; 
  GidPostProcessModes_domGAUSSnew(MESH.COOR,MESH.CN,MESH.TypeElement,BasisSTRESS,DATA.MESH.posgp,NameFileMesh,NameFile_res,DATA);

% GidPostProcessModesDOML(MESH.COOR,MESH.CN,MESH.TypeElement,BasisREACTF,DATA.MESH.posgp,NameFileMesh,NameFile_res,MESH.MaterialType,DATALOC) ;
 