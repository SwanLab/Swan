classdef Settings
    properties %optmizer access
        plotting = true
        showBC = true
        printing = true
        printing_physics = false
        monitoring = true
        monitoring_interval = 10
        maxiter = 2000
        constraint_case = 'EQUALITY';
        N_holes = [5 6 4];
        R_holes = 0.7;
        phase_holes = [0 pi/2 0];
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
        epsilon_isotropy_initial
        epsilon_isotropy_final
        HJiter0 = 1
    end
    
    properties
        perimeter = struct;
    end
    
    properties    %topopt access
        ptype
        case_file
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
        selectiveC_Cstar
        nconstr
    end
    
    methods
        function obj = Settings(case_file)
            run(case_file)
            obj.case_file=case_file;
            obj.filename = filename;
            obj.ptype = ptype;
            obj.method = method;
            obj.material = materialType;
            obj.initial_case = initial_case;
            obj.cost = cost;
            obj.weights = weights;
            obj.constraint = constraint;
            obj.nconstr = length(constraint);
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
                        
            if exist('constraint_case','var')
                obj.constraint_case = constraint_case;
            end
            
            if exist('plotting','var')
                obj.plotting = plotting;
            end
            if exist('showBC','var')
                obj.showBC = showBC;
            end
            if exist('printing','var')
                obj.printing = printing;
            end
            if exist('printing_physics','var')
                obj.printing_physics = printing_physics;
            end
            if ~obj.printing && obj.printing_physics
                warning('Physical variables will not be printed.')
            end
            if exist('monitoring','var')
                obj.monitoring = monitoring;
            end
            if exist('monitoring_interval','var')
                obj.monitoring_interval = monitoring_interval;
            end
            if exist('maxiter','var')
                obj.maxiter = maxiter;
            end
            
            if ~contains(case_file,'test','IgnoreCase',true)
                fprintf('Loaded %s: \n -Optimizer: %s \n -Cost: ',case_file,obj.optimizer)
                fprintf('%s, ',obj.cost{:})
                fprintf('\n -Constraint: ')
                fprintf('%s, ', obj.constraint{:})
                fprintf('\n -Incremental Steps: %f \n ',obj.nsteps)
            end
            
            if exist('Vfrac_final','var')
                obj.Vfrac_final = Vfrac_final;
                if ~contains(case_file,'test','IgnoreCase',true)
                    fprintf('-Volume target: %f \n ',obj.Vfrac_final)
                end
            end
            if exist('Perimeter_target','var')
                obj.Perimeter_target = Perimeter_target;
                if ~contains(case_file,'test','IgnoreCase',true)
                    fprintf('-Perimeter target: %f \n',obj.Perimeter_target)
                end
            end
            if exist('epsilon_initial','var')
                obj.epsilon_initial = epsilon_initial;
                obj.epsilon_final = epsilon_final;
            end
            if exist('HJiter0','var')
                obj.HJiter0 = HJiter0;
            end
            if exist('N_holes','var')
                obj.N_holes = N_holes;
            end
            if exist('R_holes','var')
                obj.R_holes = R_holes;
            end
            if exist('phase_holes','var')
                obj.phase_holes = phase_holes;
            end
            if exist('micro','var')
                obj.micro.alpha = micro.alpha;
                obj.micro.beta = micro.beta;
            end
            if exist('epsilon_isotropy_initial','var')
                obj.epsilon_isotropy_initial = epsilon_isotropy_initial;
                obj.epsilon_isotropy_final = epsilon_isotropy_final;
            end
            if exist('selectiveC_Cstar','var')
                obj.selectiveC_Cstar = selectiveC_Cstar;
            end
            if ~contains(filename,'test','IgnoreCase',true)
                fprintf('\n')
            end
        end
    end
end
