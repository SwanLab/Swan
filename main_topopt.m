%% path adding
addpath('./OOP FEM/')
addpath('./Topology Optimization/')

%% main
clc 
clear variables
settings=struct;
settings.method='SIMPALL';
settings.material='ISOTROPIC';
settings.topoptproblem='Compliance_st_VolPer';
HSbounds.gamma_plus=1;
HSbounds.gamma_minus=0;
HSbounds.E_plus=1;
HSbounds.E_minus=1e-3;
HSbounds.nu_plus=1/3;
HSbounds.nu_minus=1/3;

physProblem=PhysProblem_SolidMechanics('TOPOPT_TEST');
physProblem.preProcess;
rho=ones(physProblem.element.nelem,physProblem.element.geometry.ngauss);

switch settings.topoptproblem
    case 'Compliance_st_VolPer'
        test=TopOpt_Problem_Compliance_st_VolPer(rho,HSbounds,physProblem,settings);
    otherwise
        disp('Problem not added')
end

test.computeVariables;