classdef Optimizer_IPOPT < Optimizer
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'IPOPT'
    end    
        
    properties (Access = private)
        nconstr
        info
        max_iter
        constraint_tolerance
        optimality_tolerance
        cost_copy
        constraint_copy
        data
    end
    methods
        function obj = Optimizer_IPOPT(cParams)
            obj.init(cParams);
            obj.nconstr  = cParams.nconstr;
            obj.max_iter = cParams.maxiter;
            obj.niter    = -1;
        end
        function ot = get.optimality_tolerance(obj)
            ot = obj.targetParameters.optimality_tol;
        end
        
        function ct = get.constraint_tolerance(obj)
            ct = obj.targetParameters.constr_tol*1e-1;
        end
        
        function solveProblem(obj)           
            %obj.createPostProcess();
            x0 = obj.designVariable.value;
            obj.cost.computeCostAndGradient()
            funcs.objective = @(x) obj.objective(x);
            funcs.gradient = @(x) obj.gradient(x);
            funcs.constraints = @(x) obj.constraintFunction(x);
            funcs.jacobian = @(x) sparse(obj.constraint_gradient(x)');
            n = length(x0);
            funcs.jacobianstructure = @() sparse(ones(obj.nconstr,n));
            funcs.iterfunc = @(iter,fval,data) obj.outputfun_ipopt(data);
            
            options.ipopt.print_level= 0;
            options.ipopt.hessian_approximation = 'limited-memory';
            options.ipopt.limited_memory_update_type = 'bfgs';
            options.ub = ones(length(x0),1);
            options.lb = zeros(length(x0),1);
            if strcmp(obj.constraintCase,'EQUALITY')
                options.cl = zeros(obj.nconstr,1);
                options.constraintCase = 'equality';
            else
                options.cl = -Inf*ones(obj.nconstr,1); % lower bound constraint
            end
            
            options.cu = zeros(obj.nconstr,1); % upper bound constraint value
            options.ipopt.max_iter = obj.max_iter;
            options.ipopt.constr_viol_tol = obj.constraint_tolerance;
            %         options.ipopt.dual_inf_tol = optimality_tolerance;
            options.ipopt.compl_inf_tol = obj.constraint_tolerance;
            options.ipopt.tol = obj.optimality_tolerance;
            
            [x, obj.info] = ipopt(x0,funcs,options);
            
            obj.designVariable.update(x);
        end
        
    end
    
    methods (Access = private)
        
        function f = objective(obj,x)
            obj.designVariable.value = x;
            obj.cost.computeCostAndGradient()
            obj.cost_copy = obj.cost;
            f = obj.cost.value;
        end
        
        function f = constraintFunction(obj,x)
            obj.designVariable.value = x;            
            obj.constraint.computeCostAndGradient()
            obj.constraint_copy = obj.constraint;
            f = obj.constraint.value;
        end
        
        function g = gradient(obj,x)
            obj.designVariable.value = x;            
            obj.cost.computeCostAndGradient()
            obj.cost_copy=obj.cost;
            g = obj.cost.gradient;
        end
        
        function g = constraint_gradient(obj,x)
            obj.designVariable.value = x;            
            obj.constraint.computeCostAndGradient()
            obj.constraint_copy = obj.constraint;
            g = obj.constraint.gradient;
        end
        
        function stop = outputfun_ipopt(obj,data)
            stop = true;
            obj.historicalVariables.inf_du = data.inf_du;
            obj.data=data;
            obj.niter=obj.niter+1;
            obj.designVariable.update(data.x);
            obj.updateStatus();
            obj.printOptimizerVariable();
            obj.constraint_copy.lambda=zeros(obj.constraint_copy.nSF,1);
            obj.convergenceVars.reset();
            obj.convergenceVars.append(data.inf_du);
            obj.refreshMonitoring();
            obj.printHistory();
        end
        
    end
    
end