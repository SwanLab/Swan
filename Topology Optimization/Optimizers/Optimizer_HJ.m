classdef Optimizer_HJ < Optimizer_Unconstrained
    
    properties
        optimality_tol
        e2
        mean_cell_size
        % !! Move to ShFunc_Velocity (?) eventually !!
        filter
    end
    
    methods
        function obj = Optimizer_HJ(settings,epsilon,mean_cell_size)
            obj@Optimizer_Unconstrained(settings,epsilon);
            obj.e2 = settings.e2;
            obj.mean_cell_size = mean_cell_size;
            obj.max_constr_change = +Inf;
            obj.nconstr = settings.nconstr;
            % !! Move to ShFunc_Velocity (?) eventually !!
            if strcmp(settings.filter,'P1')
                settings.filter = 'PDE';
                disp('Filter P1 changed to PDE for HJ velocity regularization');
            end
            obj.filter = Filter_Boundary.create(settings);
            obj.filter.setupFromGiDFile(settings.filename,settings.ptype);
            obj.filter.preProcess;
            obj.filter.updateEpsilon(epsilon);
        end
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = (0.0175/1e-3)*obj.target_parameters.optimality_tol;
        end
        
        function phi = computeX(obj,phi,gradient)
            V = -obj.filter.regularize(phi,gradient);
            
            dt = 0.5*obj.e2*obj.line_search.kappa*obj.mean_cell_size/max(abs(V(:))) ;
            phi = obj.solvelvlset(phi,V,dt);
            
            obj.opt_cond = obj.line_search.kappa;
        end
        
        function solved_phi = solvelvlset(obj,phi,V,dt)
            for i = 1:obj.line_search.HJiter                
                phi = phi - dt*V;
            end
            solved_phi = phi;
        end
    end
end
