classdef Settings < handle
    properties %optmizer access
        plotting=true
        printing=false
        monitoring=true
        monitoring_interval=10
        maxiter=5000
    end
    properties %target parameters
        Vfrac_initial = 1;
        optimality_initial = 1e-3;
        constr_initial = 1e-3;
        Vfrac_final
        optimality_final
        constr_final
        epsilon_initial
        epsilon_final
        Perimeter_target
    end
    properties 
        perimeter=struct;
    end
    properties    %topopt access
        ptype = 'MACRO';
        filename = 'CantileverBeam_Triangle_Linear_Fine';
        method = 'SIMPALL';
        material = 'ISOTROPIC';
        initial_case = 'full';
        cost = {'compliance'}
        weights=[];
        constraint = {'volume'};
        optimizer = 'SLERP';
        kappaMultiplier
        filter = 'P1';
        TOL=struct;
        target_parameters=struct;        
        nsteps = 1; 
        micro=struct;  
        nconstr
    end
    methods
        function obj=Settings(case_file,mesh_file)
                obj.filename=mesh_file;
                
                run(case_file)               
                obj.ptype = ptype;                
                obj.method = method;
                obj.material=materialType;
                obj.initial_case = initial_case;
                obj.cost = cost;
                obj.weights=weights;
                obj.constraint = constraint;
                obj.optimizer = optimizer;
                obj.kappaMultiplier=kappaMultiplier;
                obj.filter = filterType;
                obj.nsteps = nsteps;
                if exist('TOL','var')
                    obj.TOL.rho_plus = TOL.rho_plus;
                    obj.TOL.rho_minus = TOL.rho_minus;
                    obj.TOL.E_plus = TOL.E_plus;
                    obj.TOL.E_minus = TOL.E_minus;
                    obj.TOL.nu_plus = TOL.nu_plus;
                    obj.TOL.nu_minus = TOL.nu_minus;
                end
                obj.Vfrac_initial = Vfrac_initial;
                obj.optimality_initial =optimality_initial;
                obj.constr_initial = constr_initial;
                obj.Vfrac_final = Vfrac_final;
                obj.optimality_final = optimality_final;
                obj.constr_final = constr_final;
                obj.Perimeter_target=Perimeter_target;

                
                obj.micro.alpha = micro.alpha;
                obj.micro.beta = micro.beta;
        end
    end
end
