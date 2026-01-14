function PlotModesSE_SubdomainLevel(DATA,PsiRBf,PsiDEFf,MESH)


NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
if  ~exist(NAME_MODES_FOLDER)
    mkdir(NAME_MODES_FOLDER)
end

DATA = DefaultField(DATA,'LEGEND_MODES_SE','') ; 

NAME_MODES_DISP = [NAME_MODES_FOLDER,'REACTmodes',DATA.LEGEND_MODES_SE] ;

NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [] ;
ncols = size(PsiDEFf,2) ; 
ndofs = prod(size(MESH.COOR)) ; 
BasisREACTF = zeros(ndofs,ncols) ; 
f = MESH.faceDOFSall ; 
BasisREACTF(f,:) = PsiDEFf ; 
if ~isempty(PsiRBf)
ncols = size(PsiRBf,2) ; 
BasisREACTF2 = zeros(ndofs,ncols) ; 
BasisREACTF2(f,:) = PsiRBf ; 
else
    BasisREACTF2 = [] ; 
end

BasisREACTF = [BasisREACTF2,BasisREACTF] ; 
DATALOC = []; 
GidPostProcessModesDOML(MESH.COOR,MESH.CN,MESH.TypeElement,BasisREACTF,DATA.MESH.posgp,NameFileMesh,NameFile_res,MESH.MaterialType,DATALOC) ;
 