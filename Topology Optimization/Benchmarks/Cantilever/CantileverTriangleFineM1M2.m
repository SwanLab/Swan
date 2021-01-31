filename='CantileverSquare';
%filename = 'BridgeArch';
%filename = 'Cantilever_triangle_fine';
ptype = 'MACRO';
initial_case = 'full';
cost = {'compliance'};
weights = 1;
constraint = {'volumeConstraint'};
filterType = 'P1';

Vfrac_final = 0.3;
optimality_final = 1e-3;
constr_final = 1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;

printing = false;
monitoring_interval = 1;
maxiter = 50;
printing = false;           


option = 'A';
switch option
    case 'A0'
        Vfrac_initial = 1;        
        optimizer = 'DualNestedInPrimal';
        optimizerUnconstrained = 'PROJECTED GRADIENT';
        designVariable = 'MicroParams';
        ub = 0.989;
        lb = 0.011;
        rate = 0.5;%1.05;
        nsteps = 50;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'SmoothRectangle';
        maxiter = 250;          
    case 'A'
        Vfrac_initial = 0.3;        
        optimizer = 'DualNestedInPrimal';
        optimizerUnconstrained = 'PROJECTED GRADIENT';
        designVariable = 'MicroParams';
        ub = 0.989;
        lb = 0.011;
        rate = 0.5;%1.05;
        nsteps = 1;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'Rectangle';%'SmoothRectangle';
        maxiter = 500;  
        filterType = 'P1';
        epsilon_initial = 0.01;
        epsilon_final = 0.01;
        optimality_final = 1e-5;
    case 'B'
        Vfrac_initial = 0.3;                
        optimizer = 'IPOPT';
        designVariable = 'Density';
        nsteps = 50;
        TOL.rho_plus = 1;
        TOL.rho_minus = 0;
        TOL.E_plus = 1;
        TOL.E_minus = 1e-3;
        TOL.nu_plus = 1/3;
        TOL.nu_minus = 1/3;
        method = 'SIMPALL';
        materialType = 'ISOTROPIC'; 
    case 'B2'
        Vfrac_initial = 0.3;                
        optimizer = 'MMA';
        designVariable = 'Density';
        nsteps = 1;
        TOL.rho_plus = 1;
        TOL.rho_minus = 0;
        TOL.E_plus = 1;
        TOL.E_minus = 1e-3;
        TOL.nu_plus = 1/3;
        TOL.nu_minus = 1/3;
        method = 'SIMPALL';
        materialType = 'ISOTROPIC';         
    case 'C'
        Vfrac_initial = 1;                        
        optimizer = 'IPOPT';
        designVariable = 'MicroParams';
        ub = 0.989;
        lb = 0.011;
        %rate = 1/1.05;
        nsteps = 1;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'SmoothRectangle';  
        maxiter = 30; 
    case 'D'
        Vfrac_initial = 0.3;                        
        optimizer = 'MMA';
        designVariable = 'MicroParams';
        ub = 0.989;
        lb = 0.011;
        %rate = 1/1.05;
        nsteps = 1;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'SmoothRectangle';  
        maxiter = 1000;        
        
    case 'E'
        optimizer = 'AlternatingPrimalDual';
        optimizerUnconstrained = 'PROJECTED GRADIENT';        
        designVariable = 'MicroParams';
        ub = 0.989;
        lb = 0.011;
        %rate = 1/1.05;
        nsteps = 1;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'SmoothRectangle';  
        maxiter = 1000;        
    case 'E2'
        Vfrac_initial = 1;                                
        optimizer = 'AlternatingPrimalDual';
        optimizerUnconstrained = 'PROJECTED GRADIENT';        
        designVariable = 'Density';
        nsteps = 1;
        TOL.rho_plus = 1;
        TOL.rho_minus = 0;
        TOL.E_plus = 1;
        TOL.E_minus = 1e-3;
        TOL.nu_plus = 1/3;
        TOL.nu_minus = 1/3;
        method = 'SIMPALL';
        materialType = 'ISOTROPIC'; 
                maxiter = 1000;        

            case 'E3'
        optimizer = 'DualNestedInPrimal';
        optimizerUnconstrained = 'SLERP';        
        designVariable = 'LevelSet';
        nsteps = 10;
        TOL.rho_plus = 1;
        TOL.rho_minus = 0;
        TOL.E_plus = 1;
        TOL.E_minus = 1e-3;
        TOL.nu_plus = 1/3;
        TOL.nu_minus = 1/3;
        method = 'SIMPALL';
        materialType = 'ISOTROPIC'; 
        maxiter = 200;     
end


