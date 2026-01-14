function PlotOrthogonalFunctToIntegrateCECM(Ubasis,MESH,DATA_ECM) ; 
if nargin == 0
    load('tmp1.mat')
end
% DATA_ECM = DefaultField(DATA_ECM,'PLOT_ORTHOGONAL_FUNCTIONS_TO_BE_INTEGRATED',1) ;
if DATA_ECM.PLOT_ORTHOGONAL_FUNCTIONS_TO_BE_INTEGRATED == 1
   % DATA_ECM.LabelPlot_integrand = 'BODYforces' ; 
DATA_ECM = DefaultField(DATA_ECM,'LabelPlot_integrand','IntWork') ; 
   
   
    NameFileMesh = [cd,filesep,'MODES/','Modes',DATA_ECM.LabelPlot_integrand,'.msh'] ;  
    DATALOC = [] ; 
    GidPostProcessModesFINT_LOCg(MESH.COOR,MESH.CN,MESH.TypeElement,Ubasis,MESH.posgp,NameFileMesh,DATALOC);
    
    
end