function     PlotModesGID_Interf(COORbnd,CNbREN,MESH,V,DATA,NameLoc) ;



NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;


for i  =1:size(V,2)
    V(:,i) = V(:,i)/norm(V(:,i)) ;
end


DATALOCaaa = [];
GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,V,DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOCaaa) ;