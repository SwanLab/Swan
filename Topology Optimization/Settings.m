classdef Settings < handle
    properties %optmizer access
        plotting = true
        printing = false
        monitoring = true
        monitoring_interval = 10
        maxiter = 5000
    end
    properties %target parameters
        Vfrac_initial
        optimality_initial 
        constr_initial
        Vfrac_final
        optimality_final
        constr_final
        epsilon_initial
        epsilon_final
        Perimeter_target
    end
    properties 
        perimeter = struct;
    end
    properties    %topopt access
        ptype
        filename 
        method
        material 
        initial_case 
        cost 
        weights
        constraint
        optimizer 
        kappaMultiplier
        filter 
        TOL = struct;
        target_parameters = struct;        
        nsteps
        micro = struct;  
        nconstr
    end
    methods
        function obj = Settings(case_file)       
                run(case_file)
                obj.filename = filename;
                obj.ptype = ptype;                
                obj.method = method;
                obj.material = materialType;
                obj.initial_case = initial_case;
                obj.cost = cost;
                obj.weights = weights;
                obj.constraint = constraint;
                %% RELEASE WHEN TEST --> AS BENCHMARK CASES
                % obj.nconstr = length(constraint);
                obj.optimizer = optimizer;
                obj.kappaMultiplier = kappaMultiplier;
                obj.filter = filterType;
                obj.nsteps = nsteps;                
                obj.TOL.rho_plus = TOL.rho_plus;
                obj.TOL.rho_minus = TOL.rho_minus;
                obj.TOL.E_plus = TOL.E_plus;
                obj.TOL.E_minus = TOL.E_minus;
                obj.TOL.nu_plus = TOL.nu_plus;
                obj.TOL.nu_minus = TOL.nu_minus;                
                obj.Vfrac_initial = Vfrac_initial;
                obj.optimality_initial  = optimality_initial;
                obj.constr_initial = constr_initial;
                obj.optimality_final = optimality_final;
                obj.constr_final = constr_final;
                fprintf('Loaded %s: \n -Optimizer: %s \n -Cost: %s \n -Constraint: %s \n -Incremental Steps: %f \n ',...
                    case_file,obj.optimizer,char(obj.cost),char(obj.constraint),obj.nsteps)
                if exist('Vfrac_final','var')
                    obj.Vfrac_final = Vfrac_final;
                    fprintf('-Volume target: %f \n ',obj.Vfrac_final)
                end
                if exist('Perimeter_target','var')
                    obj.Perimeter_target = Perimeter_target;
                    fprintf('-Perimeter target: %f \n',obj.Perimeter_target)
                end
                if exist('micro','var')
                    obj.micro.alpha = micro.alpha;
                    obj.micro.beta = micro.beta;
                end
                fprintf('\n')
        end
    end
end
