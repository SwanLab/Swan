function [matCh,h_C_0] = evaluate_laminates(phi_date,printi,physical_type,irun)

Vfrac = 0.60;

firstmesh = 1;
penalty = 0.1;
lambda0 = 0;


lagrange_type = 'AUGMENTED';
function_type = 'MIN_STIFF_INV';


%snaps = [6];


%file_path = ['/home/aferrer/Documents/Doctorat/Tesi/MicroEscalaDT/Vademecum/',num2str(irun),'/'];
file_path = [pwd,'/VademecumThermal_LaminadoNoPeriodico/',num2str(irun),'/'];
file_write = file_path;
inicial = 1;
%unix(['rm -r ',file_path]);
if inicial
    unix(['mkdir ',file_path]);
end
phi_v = phi_date;
theta_v = 0;


% data file name, read data 
TYPE = 'MICRO'; %'MACRO' 'MICRO'

%path = pwd;
%addpath(genpath('/home/aferrer/Documents/Doctorat/Tesi/VortexForcesAirfoil'))
%AirfoilGidData([path,'/Results/RVE04N3'],0);
addpath(genpath('/home/aferrer/Documents/Doctorat/Tesi/ElasticThermal_Regularized'))
path_mesh = pwd;
file_name = 'RVE04N3';%'HeatTransferMacro';%'RVE04N3';%''Cantiliberbeam'%'BridgeSimetric';%'Bridge';%'Airfoil_min_conditions';%'Cantiliberbeam_antiHolesMiddle';%'Cantiliberbeam_antiHolesMiddle';%'CantiliberbeamHolesMiddle';%'Cantiliberbeam';%'Cantiliberbeam_anti';%'Cantiliberbeam';%'BendingBeam';%'Airfoil_min_conditions';%'Cantiliberbeam';%'Airfoil_min_conditions';%'Airfoil_skin_scan_rib_with_holes';;%'Cantiliberbeam';%'Airfoil_skin_scan_rib';%'Cantiliberbeam';%'Airfoil';%'Cantiliberbeam';%'Airfoil_currado'; %'Cantiliberbeam'
path_mesh = [path_mesh,'/Results/RVE04N3/',file_name,'.gid'];
unix(['cp ',path_mesh,'/',file_name,'.m ',path_mesh,'/',file_name,'_MESH_1.m']); 

%file_name = 'CANTILEVER_TRI';%% Perfil_Alar %CANTILEVER_TRI %CANTILEVER_TRI_128 %_25600';  CANTILEVER_TRI CANTILEVER_QUAD % 'RVE04N3'; 
file_gid = [file_path,file_name];%_25600';   % 'RVE04N3'; 

% MACRO STRAIN
switch physical_type
    case 'ELASTIC'
        alpha = phitheta2strain(phi_v,theta_v)';
    case 'THERMAL'
        alpha = [cos(phi_v) sin(phi_v)]';
end
weight = alpha*alpha'/(alpha'*alpha);

%First mesh
imesh = firstmesh;
[fext,element,fixnodes,problembsc,coordinates,phifunct_n,dim,Msmooth,Stiff_smooth,emass,Group] = call_read_data(file_name,Vfrac,imesh,lagrange_type,function_type,file_write,TYPE,physical_type,phi_v,theta_v);
%fext.pointload(:,3) = fext.pointload(:,3)/100;

if problembsc.flag_change_micro == 1
element = vademecum_info(element,dim,phyisical_type);
else
element.Ch = 0;
end
%Group = create_group(coordinates,dim,element,4);
if ~isempty(Group)
    element.Group = repmat(Group(:,2),1,element.ngaus)';
else
    Group = [[1:dim.nelem]',zeros(dim.nelem,1)];
    element.Group = repmat(Group(:,2),1,element.ngaus)';
end


% Initial value of 
h_C_0 = 1;
Perim0 = 1;
%h_C_0 = 8.527093766160089;

u0 = zeros(dim.nndof,1);
lambda_n = lambda0;
penalty_n = penalty;
flag_up_date_lambda = 1;
flag_change_micro_ini = 0;
[cost_n,theta_n,constr_n,vol_n,h_C_0,~,~,structural_values,gfunc_til_n,g_nodal_n,norm_g_n,g_ortho,norm_g_ortho,g_phi,d_u,vdisp,nbdata,post,phi_gaus,theta_gaus,~,~,Perim0] = ...
equilibrium_update(phifunct_n,weight,h_C_0,Perim0,lambda_n,penalty_n,element,problembsc,fixnodes,coordinates,fext,dim,Msmooth,Stiff_smooth,emass,0,flag_change_micro_ini,Group,0,file_gid,file_write,u0);
cost_n = h_C_0/h_C_0 +lambda_n*constr_n + 0.5*1/penalty_n*constr_n^2+ problembsc.alpha_perimeter*Perim0/Perim0;

element.phi_gaus = phi_gaus;
element.theta_gaus = theta_gaus;

if isempty(fext.pointload)
    structural_values.fext = [];
else
structural_values.fext = change_Fext_format(fext.pointload);
end

if printi 
printinfo(1,problembsc,d_u,file_gid,0,coordinates,element,dim,phifunct_n,g_nodal_n,gfunc_til_n,structural_values,post,g_ortho,norm_g_n,vdisp,nbdata);
  % postproceso y almacenaje
print_variables2gid_and_screen(0,cost_n,theta_n,1,vol_n,lambda_n,penalty_n,problembsc,d_u,file_gid,file_write,coordinates,element,...
        dim,phifunct_n,g_nodal_n/g_phi,gfunc_til_n,post,g_ortho,vdisp,nbdata,imesh,h_C_0/h_C_0,constr_n,structural_values,Perim0/Perim0,norm_g_n); 
end

matCh = structural_values.matCh;
