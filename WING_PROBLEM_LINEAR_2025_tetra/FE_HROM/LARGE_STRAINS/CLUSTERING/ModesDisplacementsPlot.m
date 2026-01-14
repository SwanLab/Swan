function ModesDisplacementsPlot(DATAoffline,NAME_BASE,MESH,BasisU) 

% PLOTTING DISPLACEMENT MODES
% ***************************
if DATAoffline.PLOT_MODES_WITH_DOFR == 1
    
    NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
    if  ~exist(NAME_MODES_FOLDER)
        mkdir(NAME_MODES_FOLDER)
    end
    NAME_MODES_DISP = [NAME_MODES_FOLDER,NAME_BASE,'DISPmodes'] ;
    
    NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
    NameFile_res = [NAME_MODES_DISP,'.res'] ;
    DATALOC = [] ;
    GidPostProcessModesDOML(MESH.COOR,MESH.CN,MESH.TypeElement,BasisU,DATA.MESH.posgp,NameFileMesh,NameFile_res,MESH.MaterialType,DATALOC) ;
    
    
    disp('Set PLOT_MODES_WITH_DOFR = 0 to generate the modes. This is just for post-process purposes' )
    return
end