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
        
        function obj = Optimizer_AugLag(settings)
            
            ocS.settings       = settings;
            ocS.designVariable = settings.designVar;
            ocS.monitoring     = settings.monitoring;
            obj@Optimizer_Constrained(ocS);
            
            
            augLagS.nconstr = settings.nconstr;
            augLagS.constraintCase = settings.constraint_case;
            
            obj.augLagrangian = AugmentedLagrangian(augLagS);
            obj.optimizer_unconstr = settings.optimizer_unconstr;
            
            obj.lambda  = zeros(1,settings.nconstr);
            obj.penalty = ones(1,settings.nconstr);
            
        end
        
        function x = update(obj,x0,cost,constraint)
            obj.augLagrangian.link(cost,constraint);
            obj.optimizer_unconstr.target_parameters = obj.target_parameters;
            
            obj.updateDualVariable(constraint);
            obj.augLagrangian.updateBecauseOfDual(obj.lambda,obj.penalty);
            
            obj.updatePrimalVariable(x0);
            
            obj.updateConvergenceStatus(constraint);
            
            x = obj.x;
        end
        
    end
    
    methods (Access = private)
        
        function updatePrimalVariable(obj,x0)
            obj.optimizer_unconstr.init(x0,obj.augLagrangian);
            while ~obj.optimizer_unconstr.has_converged
                x = obj.optimizer_unconstr.update(x0,obj.augLagrangian);
                obj.stop_vars = obj.optimizer_unconstr.stop_vars;
            end
            
            if ~obj.optimizer_unconstr.designImproved
                x = x0;
            end
            
            obj.x = x;
        end
    
        function updateConvergenceStatus(obj,constraint)
            active_constr = obj.penalty > 0;
            hasNotConverged = obj.optimizer_unconstr.opt_cond >=  obj.optimizer_unconstr.optimality_tol || any(any(abs(constraint.value(active_constr)) > obj.optimizer_unconstr.constr_tol(active_constr)));
            obj.has_converged = ~hasNotConverged;
        end
        
        function updateDualVariable(obj,constraint)
            l   = obj.lambda;
            rho = obj.penalty;
            c   = constraint.value';
            l = l + rho.*c;
            obj.lambda = l;
        end
        
    end
    
end