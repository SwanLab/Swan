filename='Cantilever_triangle_coarse';
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

option = 'D';
switch option
    case 'A'
        optimizer = 'DualNestedInPrimal';
        optimizerUnconstrained = 'PROJECTED GRADIENT';
        designVariable = 'MicroParams';
        ub = 0.98;
        lb = 0.02;
        kfrac = 1.05;
        nsteps = 50;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'SmoothRectangle';
    case 'B'
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
    case 'C'
        optimizer = 'IPOPT';
        designVariable = 'MicroParams';
        ub = 0.98;
        lb = 0.02;
        %kfrac = 1.05;
        nsteps = 20;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'SmoothRectangle';  
        maxiter = 10;        
    case 'D'
        optimizer = 'MMA';
        designVariable = 'MicroParams';
        ub = 0.98;
        lb = 0.02;
        %kfrac = 1.05;
        nsteps = 20;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'SmoothRectangle';  
        maxiter = 100;        
        
    case 'E'
        optimizer = 'AlternatingPrimalDual';
        optimizerUnconstrained = 'PROJECTED GRADIENT';        
        designVariable = 'MicroParams';
        ub = 0.98;
        lb = 0.02;
        %kfrac = 1.05;
        nsteps = 20;
        homegenizedVariablesComputer = 'ByVademecum';
        vademecumFileName = 'SmoothRectangle';  
        maxiter = 10;        
end


