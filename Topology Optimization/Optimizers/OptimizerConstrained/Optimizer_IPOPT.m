classdef Optimizer_IPOPT < Optimizer
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'IPOPT'
    end
    
    properties (Access = private)
        nconstr
        info
        max_iter
        constraintTolerance
        optimalityTolerance
        data
        upperBound
        lowerBound
        functions
        nX
        options
    end
    
    methods (Access = public)
        
        function obj = Optimizer_IPOPT(cParams)
            obj.init(cParams);
            obj.upperBound = cParams.uncOptimizerSettings.ub;
            obj.lowerBound = cParams.uncOptimizerSettings.lb;
            obj.nconstr    = cParams.nconstr;
            obj.max_iter   = cParams.maxiter;
            obj.niter      = -1;
            obj.nX         = length(obj.designVariable.value);
            obj.createFunctions();
            obj.createOptions();
        end
        
        function solveProblem(obj)
            obj.cost.computeCostAndGradient();
            obj.updateIpoptOptions();
            x = obj.callIpopt();
            obj.designVariable.update(x);
        end
        
    end
    
    methods (Access = private)
        
        function x = callIpopt(obj)
            x0  = obj.designVariable.value;
            fun = obj.functions;
            opt = obj.options;
            [x, obj.info] = ipopt(x0,fun,opt);
        end
        
        function createFunctions(obj)
            funcs.objective         = @(x) obj.objective(x);
            funcs.gradient          = @(x) obj.gradient(x);
            funcs.constraints       = @(x) obj.constraintFunction(x);
            funcs.jacobian          = @(x) sparse(obj.constraint_gradient(x)');
            funcs.jacobianstructure = @() sparse(ones(obj.nconstr,obj.nX));
            funcs.iterfunc          = @(iter,fval,data) obj.outputfun_ipopt(data);
            obj.functions           = funcs;
        end
        
        function f = objective(obj,x)
            obj.designVariable.value = x;
            obj.cost.computeCostAndGradient()
            f = obj.cost.value;
        end
        
        function f = constraintFunction(obj,x)
            obj.designVariable.value = x;
            obj.constraint.computeCostAndGradient()
            f = obj.constraint.value;
        end
        
        function g = gradient(obj,x)
            obj.designVariable.value = x;
            obj.cost.computeCostAndGradient()
            g = obj.cost.gradient;
        end
        
        function g = constraint_gradient(obj,x)
            obj.designVariable.value = x;
            obj.constraint.computeCostAndGradient()
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
            obj.dualVariable.value = zeros(obj.constraint.nSF,1);
            
            obj.convergenceVars.reset();
            obj.convergenceVars.append(data.inf_du);
            obj.refreshMonitoring();
            obj.printHistory();
        end
        
        function createOptions(obj)
            opt.ipopt.print_level= 0;
            opt.ipopt.hessian_approximation = 'limited-memory';
            opt.ipopt.limited_memory_update_type = 'bfgs';
            opt.ub = obj.upperBound*ones(obj.nX,1);
            opt.lb = obj.lowerBound*ones(obj.nX,1);
            if strcmp(obj.constraintCase,'EQUALITY')
                opt.cl = zeros(obj.nconstr,1);
                opt.constraintCase = 'equality';
            else
                opt.cl = -Inf*ones(obj.nconstr,1);
            end
            opt.cu                    = zeros(obj.nconstr,1);
            opt.ipopt.max_iter        = obj.max_iter;
            opt.ipopt.constr_viol_tol = obj.obtainConstraintTolerance();
            opt.ipopt.compl_inf_tol   = obj.obtainConstraintTolerance();
            opt.ipopt.tol             = obj.obtainOptimalityTolerance();
            obj.options = opt;
        end
        
        function updateIpoptOptions(obj)
            obj.ipopt.constr_viol_tol = obj.obtainConstraintTolerance();
            obj.ipopt.compl_inf_tol   = obj.obtainConstraintTolerance();
            obj.ipopt.tol             = obj.obtainOptimalityTolerance();            
        end
        
        function ot = obtainOptimalityTolerance(obj)
            ot = obj.targetParameters.optimality_tol;
        end
        
        function ct = obtainConstraintTolerance(obj)
            ct = obj.targetParameters.constr_tol*1e-1;
        end
        
    end
    
end