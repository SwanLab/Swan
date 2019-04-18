classdef Optimizer_AugLag < Optimizer_Constrained
    
    properties (GetAccess = public, SetAccess = private)
        optimizer_unconstr
        augLagrangian
        penalty
    end
    
    properties (Access = private)
        x
        lambda
    end
    
    methods (Access = public)
        
        function obj = Optimizer_AugLag(cParams)            
            obj.init(cParams);
            
            obj.convergenceVars = cParams.optimizer_unconstr.convergenceVars;
            augLagS.constraintCase = obj.constraintCase;
            
            obj.augLagrangian = AugmentedLagrangian(augLagS);
            obj.augLagrangian.link(obj.cost,obj.constraint);
            obj.optimizer_unconstr = cParams.optimizer_unconstr;
            
            nConstraints = obj.constraint.nSF;
            obj.lambda  = zeros(1,nConstraints);
            obj.penalty = ones(1,nConstraints);
        end
        
        function x = update(obj,x0)
            obj.optimizer_unconstr.target_parameters = obj.target_parameters;            
            obj.updateDualVariable();
            obj.augLagrangian.updateBecauseOfDual(obj.lambda,obj.penalty);            
            obj.updatePrimalVariable(x0);            
            obj.updateConvergenceStatus();            
            x = obj.x;
        end
        
    end
    
    methods (Access = private)
        
        function updatePrimalVariable(obj,x0)
            obj.optimizer_unconstr.init(x0,obj.augLagrangian);
            while ~obj.optimizer_unconstr.hasConverged
                x = obj.optimizer_unconstr.update(x0,obj.augLagrangian);
            end
            
            if ~obj.optimizer_unconstr.designImproved
                x = x0;
            end
            
            obj.x = x;
        end
    
        function updateConvergenceStatus(obj)
            active_constr = obj.penalty > 0;
            isNotOptimal  = obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol;
            isNotFeasible = any(any(abs(obj.constraint.value(active_constr)) > obj.optimizer_unconstr.constr_tol(active_constr)));
            hasNotConverged = isNotOptimal || isNotFeasible;
            obj.hasConverged = ~hasNotConverged;
        end
        
        function updateDualVariable(obj)
            l   = obj.lambda;
            rho = obj.penalty;
            c   = obj.constraint.value';
            l = l + rho.*c;
            obj.lambda = l;
        end
        
    end
    
end