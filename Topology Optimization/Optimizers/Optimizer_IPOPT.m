classdef Optimizer_IPOPT < Optimizer_Constrained
    properties
        m
        info
        max_iter
        constraint_tolerance
        optimality_tolerance
        cost_copy
        constraint_copy
        data        
    end
    methods
        function obj = Optimizer_IPOPT(settings,mesh)
            obj@Optimizer_Constrained(settings,mesh,settings.monitoring);
            obj.m = settings.nconstr;
            obj.max_iter = settings.maxiter;
            obj.niter=-1;
        end
        function optimality_tolerance = get.optimality_tolerance(obj)
            optimality_tolerance = obj.target_parameters.optimality_tol;
        end
        function constraint_tolerance = get.constraint_tolerance(obj)
            constraint_tolerance = obj.target_parameters.constr_tol*1e-1;
        end
        
        function x = solveProblem(obj,x_ini,cost,constraint,istep,nstep)
            obj.createPostProcess(cost,constraint);
            cost.computeCostAndGradient(x_ini)
            funcs.objective = @(x) obj.objective(x,cost);
            funcs.gradient = @(x) obj.gradient(x,cost);
            funcs.constraints = @(x) obj.constraint(x,constraint);
            funcs.jacobian = @(x) sparse(obj.constraint_gradient(x,constraint)');
            n = length(x_ini);
            funcs.jacobianstructure = @() sparse(ones(obj.m,n));
            funcs.iterfunc = @(iter,fval,data) obj.outputfun_ipopt(data,istep,nstep);
            
            options.ipopt.print_level= 0;
            options.ipopt.hessian_approximation = 'limited-memory';
            options.ipopt.limited_memory_update_type = 'bfgs';
            options.ub = ones(length(x_ini),1);
            options.lb = zeros(length(x_ini),1);
            if strcmp(obj.constraint_case,'EQUALITY')
                options.cl = zeros(obj.m,1);
                options.constraint_case = 'equality';
            else
                options.cl = -Inf*ones(obj.m,1); % lower bound constraint
            end
            
            options.cu = zeros(obj.m,1); % upper bound constraint value
            options.ipopt.max_iter = obj.max_iter;
            options.ipopt.constr_viol_tol = obj.constraint_tolerance;
            %         options.ipopt.dual_inf_tol = optimality_tolerance;
            options.ipopt.compl_inf_tol = obj.constraint_tolerance;
            options.ipopt.tol = obj.optimality_tolerance;
            
            [x, obj.info] = ipopt(x_ini,funcs,options);
        end
        
    end
    methods 
        function f = objective(obj,x,cost)
            cost.computeCostAndGradient(x)
            obj.cost_copy=cost;
            f = cost.value;
        end
        function f = constraint(obj,x,constraint)
            constraint.computeCostAndGradient(x)
            obj.constraint_copy=constraint;
            f = constraint.value;
        end
        function g = gradient(obj,x,cost)
            cost.computeCostAndGradient(x)   
            obj.cost_copy=cost;
            g = cost.gradient;
        end
        function g = constraint_gradient(obj,x,constraint)
            constraint.computeCostAndGradient(x)    
            obj.constraint_copy=constraint;
            g = constraint.gradient;
        end
        function stop = outputfun_ipopt(obj,data,istep,nstep)            
            stop = true;
            obj.data=data;
            obj.niter=obj.niter+1;
            obj.print(data.x,obj.niter,obj.cost_copy,obj.constraint_copy);            
            obj.constraint_copy.lambda=zeros(obj.constraint_copy.nSF,1);
            obj.monitoring.refresh(data.x,obj.niter,obj.cost_copy,obj.constraint_copy,data.inf_du,obj.has_converged || obj.niter > obj.maxiter*(istep/nstep),istep,nstep);            
            obj.writeToFile(istep,obj.cost_copy,obj.constraint_copy)
        end
    end
end