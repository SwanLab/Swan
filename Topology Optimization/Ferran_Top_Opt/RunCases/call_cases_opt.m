function call_cases_opt
clear;clc;clear global;close all;
bpause = false; % 1/0 pauses after each run
bnotification = false; % 1/0 notificates after each run
Runner = char(java.net.InetAddress.getLocalHost.getHostName);
switch Runner
    case 'FERRY58-PC'
        working_path = 'E:\Dropbox\Ferran\CodeTopOpt';
        results_path = 'E:\TFM\42 - User manual';
    case 'FERRY58-PORT'
        working_path = 'D:\Dropbox\Ferran\CodeTopOpt';
        results_path = 'D:\TFM\TestingMMA';
    case 'alex'
        working_path = '/home/alex/Dropbox/Ferran/CodeTopOpt';
        results_path = '/home/alex/Desktop/Results/Results';
    case 'GOKU'
        working_path = '/home/aferrer/Dropbox/Ferran/CodeTopOpt';
        results_path = '/home/aferrer/Desktop/ResultsMicro/Results';
    otherwise
        error('Runner not detected.');
end
top_opt_path = []; addpath(genpath(working_path));
if ispc
    rmpath([working_path,'\ipopt\ipoptLinux']); 
else
    rmpath([working_path,'/ipopt/ipoptWindows']);
end
% format long
% top_opt_path = working_path; % this option deletes the path after running
% file_name = 'CantiliberbeamSymFineFineNew'    'CantiliberbeamSymNew'
%             'GrippingNew'
%             'BridgeSimetricNew'
%             'SquareNew'
%             'BicycleNew'              'BicycleFineNew'
%             'RVE_Square'              'RVE_Square_Fine'          'RVE_Square_FineFine'
%             'RVE_Hexagonal'           'RVE_Hexagonal_Fine'       'RVE_Hexagonal_FineFine' 
%             'TrencalosSupport'        'TrencalosSupportFine'
%  for more cases go to \CodeTopOpt\Meshes (some may not be compatible
%  with current version of the code)

%% METHODS
interpolation_method='SIMP_ALL'; % 'SIMP' 'SIMP_ALL' 'all'
kernel='P1'; % 'P1' 'P0' 'PDE' 'all'
algorithm='LS'; % 'PG' 'LS' 'MMA' 'fmincon' 'GCMMA' 'IPOPT' 'all' 'check_derivative'
TYPE = 'MACRO'; % 'MACRO' 'MICRO'
problem_type = 'min_compliance+lam*per_st_vol'; % ***see create_optimization_problem.m for options***
perimeter_case = 'domain_and_contour_perimeter'; % 'domain_perimeter' 'domain_and_contour_perimeter'

%% OPTIMIZATION DATA
% problem2run = read_cases('cases2run.xlsx'); % input from excel table

if ~exist('problem2run','var')
    save_directory = {results_path};
    name = {'ExampleGID'}; % can be cell of strings (without spaces) or numeric array
    mesh={'CantileverExampleProblem'};
    
    Vf=0.3;
    Ptarget=1;

    % Material data
    E_plus=1;           nu_plus=1/3;
    E_minus=1e-3;      nu_minus=1/3;

    % Generate table
    problem2run = initialize_table(name,mesh,Vf,Ptarget,E_plus,E_minus,nu_plus,nu_minus,save_directory);
end
init_design_variable = 'full'; % 'square' 'circle' 'horiz' 'rand' 'full' 'feasible'
incropt = initialize_incremental_optimization(TYPE);

run_cases_opt (interpolation_method,kernel,algorithm,TYPE,problem_type,incropt,init_design_variable,perimeter_case,problem2run,top_opt_path,bpause,bnotification);
end