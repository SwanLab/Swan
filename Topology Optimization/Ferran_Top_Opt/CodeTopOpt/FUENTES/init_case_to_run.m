function  [file_write,file_name] = init_case_to_run (file_name)

if nargin == 0
    file_name = 'CantiliberbeamSym';%'GrippingFineFine';%'CantiliberbeamSymStructuredFineFine';%GrippingFineFine';%'Gripping';%'Gripping';
                %'Gripping';
                %'Cantiliberbeam2Holes';%'CantiliberbeamSym';%'Square';%;'CantiliberbeamSymFineFine';%'CantiliberbeamSym';
                %'CantiliberbeamSym';%'CantiliberbeamSym';%'Cantiliberbeam2Holes';%'CantiliberbeamSym'; %'CantiliberbeamNoSymNoHoles';%'CantiliberbeamBandNoSym';%'CantiliberbeamBand';%'CantiliberbeamManyHoles';
                                     %'Cantiliberbeam2Holes';%'Cantiliberbeam8Holes';%'CantiliberbeamNoSymNoHoles';%'CantiliberbeamNoSymNoHolesFine';
                                     %'CantiliberbeamNoSymNoHolesFineLaminate';%'CantiliberbeamNoSymNoHolesFine';%'CantiliberbeamNoSymNoHolesFineHoriz';
                                     %'BendingBeam';%'CantiliberbeamNoSymNoHoles';%'Bridge';%'BridgeSimetric';%'RVE04N3';%'OrtotropicHeatConductor';
                                     %'HeatTransferMacro';%'RVE_Square';%''Cantiliberbeam'%'BridgeSimetric';%'Bridge';%'Airfoil_min_conditions';
                                     %'Cantiliberbeam_antiHolesMiddle';%'Cantiliberbeam_antiHolesMiddle';%'CantiliberbeamHolesMiddle';
                                     %'Cantiliberbeam';%'Cantiliberbeam_anti';%'Cantiliberbeam';%'BendingBeam';%'Airfoil_min_conditions';
                                     %'Cantiliberbeam';%'Airfoil_min_conditions';%'Airfoil_skin_scan_rib_with_holes';;%'Cantiliberbeam';
                                     %'Airfoil_skin_scan_rib';%'Cantiliberbeam';%'Airfoil';%'Cantiliberbeam';%'Airfoil_currado';
                                     %%'Cantiliberbeam';
end

% oper_sys = 'windows';                                 
%  switch oper_sys
%      case 'windows'
%          path_mesh = 'D:\TFM\Codi Opt Top';
%          path_mesh_read = [path_mesh,'\Results\RVE04N3\',file_name,'.gid'];
%          file_write = [path_mesh,'\Vademecum\','1','\',file_name];
%          copyfile([path_mesh_read,'\',file_name,'.m'],[path_mesh_read,'\',file_name,'_MESH_1.m'])
%      case 'linux'
%          path_mesh = '/home/aferrer/Dropbox/Amstutz/CodigoIntegratedLevelSetGradientPerimeter/';
%          path_mesh_read = [path_mesh,'/Results/RVE04N3/',file_name,'.gid'];
%          file_write = [path_mesh,'/Vademecum/','1','/',file_name];
%          copyfile([path_mesh_read,'/',file_name,'.m'],[path_mesh_read,'/',file_name,'_MESH_1.m'])
%  end

if ispc % WINDOWS
    path_mesh = fileparts(fileparts(mfilename('fullpath')));
    path_mesh_read = [path_mesh,'\Meshes\',file_name,'.gid'];
    file_write = [path_mesh,'\Vademecum\','1','\',file_name];
    copyfile([path_mesh_read,'\',file_name,'.m'],[path_mesh_read,'\',file_name,'_MESH_1.m'])
    
elseif isunix % LINUX
    path_mesh = fileparts(fileparts(mfilename('fullpath')));
    path_mesh_read = [path_mesh,'/Meshes/',file_name,'.gid'];
    file_write = [path_mesh,'/Vademecum/','1','/',file_name];
    copyfile([path_mesh_read,'/',file_name,'.m'],[path_mesh_read,'/',file_name,'_MESH_1.m'])
    
 %   path_mesh = '/home/aferrer/Dropbox/Amstutz/CodigoIntegratedLevelSetGradientPerimeter/';
 %   path_mesh_read = [path_mesh,'/Results/RVE04N3/',file_name,'.gid'];
 %   file_write = [path_mesh,'/Vademecum/','1','/',file_name];
 %   copyfile([path_mesh_read,'/',file_name,'.m'],[path_mesh_read,'/',file_name,'_MESH_1.m'])
    
else
    error('Operating system not included!.');
end

end