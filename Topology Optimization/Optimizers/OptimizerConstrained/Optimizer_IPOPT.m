classdef Optimizer_IPOPT < Optimizer
    
    properties (GetAccess = public, SetAccess = protected)
        type = 'IPOPT'
    end
    
    properties (Access = private)
        nConstr
        info
        constraintTolerance
        optimalityTolerance
        data
        upperBound
        lowerBound
        functions
        nX
        options
        constraintCase
    end
    
    methods (Access = public)
        
        function obj = Optimizer_IPOPT(cParams)
            obj.initOptimizer(cParams);
            obj.upperBound = cParams.uncOptimizerSettings.ub;
            obj.lowerBound = cParams.uncOptimizerSettings.lb;
            obj.nConstr    = cParams.constraint.nSF;
            obj.maxIter    = cParams.maxIter;
            obj.constraintCase = cParams.constraintCase;
            obj.nIter      = -1;
            obj.nX         = length(obj.designVariable.value);
            obj.createFunctions();
            obj.createOptions();
            obj.outputFunction.monitoring.create(cParams);            
        end
        
        function solveProblem(obj)
            obj.cost.computeFunctionAndGradient();
            obj.updateIpoptOptions();
            obj.designVariable.updateOld();
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
            funcs.jacobianstructure = @() sparse(ones(obj.nConstr,obj.nX));
            funcs.iterfunc          = @(iter,fval,data) obj.outputfun_ipopt(data);
            obj.functions           = funcs;
        end
        
        function f = objective(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            f = obj.cost.value;
        end
        
        function f = constraintFunction(obj,x)
            obj.designVariable.update(x);
            obj.constraint.computeFunctionAndGradient()
            f = obj.constraint.value;
        end
        
        function g = gradient(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient()
            g = obj.cost.gradient;
        end
        
        function g = constraint_gradient(obj,x)
            obj.designVariable.update(x);
            obj.constraint.computeFunctionAndGradient()
            g = obj.constraint.gradient;
        end
        
        function stop = outputfun_ipopt(obj,data)
             obj.nIter = obj.nIter+1;


             s.hasFinished          = false;%data.inf_pr;obj.hasFinished;
             s.nIter                = obj.nIter;
             s.inf_pr               = data.inf_pr;%obj.KKTnorm;
             s.inf_du               = data.inf_du;   
             s.outitFrac            = obj.nIter; %obj.outit/obj.maxoutit;
             obj.designVariable.updateOld();
                obj.outputFunction.monitoring.compute(s);
            
%             obj.historicalVariables.inf_pr = data.inf_pr;
%             obj.historicalVariables.inf_du = data.inf_du;
%             obj.data = data;
%             obj.nIter = obj.nIter+1;
%             if ~isempty(data.x)
%                 obj.designVariable.update(data.x);
%                 normXsquare = obj.designVariable.computeL2normIncrement();
%                 obj.designVariable.updateOld();
%                 incX = sqrt(normXsquare);
%             end
%             obj.updateStatus();
%             obj.printOptimizerVariable();
%             obj.dualVariable.value = zeros(obj.constraint.nSF,1);
%             
%             
%             obj.convergenceVars.reset();
%             obj.convergenceVars.append(data.inf_pr);
%             obj.convergenceVars.append(data.inf_du);
%             obj.convergenceVars.append(incX);
%             obj.refreshMonitoring();
%             obj.printHistory();
             stop = true;

        end
        
        function createOptions(obj)
            opt.ipopt.print_level= 0;
            opt.ipopt.hessian_approximation = 'limited-memory';
            opt.ipopt.limited_memory_update_type = 'bfgs';
            opt.ipopt.mu_strategy                = 'adaptive'; 
            opt.ipopt.mu_init                = 0.01;
            opt.ub = obj.upperBound*ones(obj.nX,1);
            opt.lb = obj.lowerBound*ones(obj.nX,1);
            if strcmp(obj.constraintCase,'EQUALITY')
                opt.cl = zeros(obj.nConstr,1);
                opt.constraintCase = 'equality';
            else
                opt.cl = -Inf*ones(obj.nConstr,1);
            end
            opt.cu                    = zeros(obj.nConstr,1);
            opt.ipopt.max_iter         = obj.maxIter;
            opt.ipopt.constr_viol_tol = obj.obtainConstraintTolerance();
            opt.ipopt.compl_inf_tol   = obj.obtainConstraintTolerance();
            opt.ipopt.tol             = obj.obtainOptimalityTolerance();
            obj.options = opt;
        end
        
        function updateIpoptOptions(obj)
            obj.options.constr_viol_tol = obj.obtainConstraintTolerance();
            obj.options.compl_inf_tol   = obj.obtainConstraintTolerance();
            obj.options.tol             = obj.obtainOptimalityTolerance();
        end
        
        function ot = obtainOptimalityTolerance(obj)
            ot = obj.targetParameters.optimality_tol;
        end
        
        function ct = obtainConstraintTolerance(obj)
            constrTol = obj.targetParameters.constr_tol; 
            ct = 0.1*constrTol;
        end
        
    end
    
end