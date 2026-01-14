function PlotModesDEF_SubdomainLevel(DATA,Phi,MESH)


NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
if  ~exist(NAME_MODES_FOLDER)
    mkdir(NAME_MODES_FOLDER);
end

 
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,'DEFmodes'] ;

NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
% DATALOC = [] ;
% ncols = size(PsiDEFf,2) ; 
% ndofs = prod(size(MESH.COOR)) ; 
% BasisREACTF = zeros(ndofs,ncols) ; 
% f = MESH.faceDOFSall ; 
% BasisREACTF(f,:) = PsiDEFf ; 
% ncols = size(PsiRBf,2) ; 
% BasisREACTF2 = zeros(ndofs,ncols) ; 
% BasisREACTF2(f,:) = PsiRBf ; 
% 
% BasisREACTF = [BasisREACTF2,BasisREACTF] ; 
DATALOC = []; 

for i = 1:size(Phi,2) 
    Phi(:,i) = Phi(:,i)/norm(Phi(:,i)) ;  
end

GidPostProcessModesDOML(MESH.COOR,MESH.CN,MESH.TypeElement,Phi,DATA.MESH.posgp,NameFileMesh,NameFile_res,MESH.MaterialType,DATALOC) ;
 