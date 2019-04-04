classdef Optimizer_Unconstrained < Optimizer
    
    properties (Access = public)
        objfunc
        opt_cond
    end
    
    properties (GetAccess = public, SetAccess = protected)
        designImproved
        max_constr_change
    end
    
    properties (GetAccess = public, SetAccess = private)
        line_search
        scalar_product
        constr_tol
    end
    
    methods (Access = public, Abstract)
        
        computeX(obj)
        
    end
    
    
    methods (Access = public)
        
        function obj = Optimizer_Unconstrained(settings,epsilon)
            obj@Optimizer(settings);
            obj.line_search = LineSearch.create(settings,epsilon);
            obj.scalar_product = ScalarProduct(settings.filename,epsilon);
        end
        
        function x = updateX(obj,x_ini,cost,constraint)
            x = obj.computeX(x_ini,obj.objfunc.gradient);
            
            obj.updateObjectiveFunction(x,cost,constraint);
            
            obj.updateConvergenceParams(x,x_ini);
            
            if ~obj.has_converged
                obj.line_search.computeKappa();
            end
        end
        
    end
    
    methods (Access = private)
        
        function updateObjectiveFunction(obj,x,cost,constraint)
            cost.computeCostAndGradient(x);
            constraint.computeCostAndGradient(x);
            constraint = obj.setConstraint_case(constraint);
            obj.objfunc.computeFunction(cost,constraint)
        end
        
        function updateConvergenceParams(obj,x,x_ini)
            incrNormL2  = obj.norm_L2(x,x_ini);
            incrCost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
            
            obj.designImproved = incrCost < 0 && incrNormL2 < obj.max_constr_change;
            
            obj.has_converged = obj.designImproved || obj.line_search.kappa <= obj.line_search.kappa_min;
            
            obj.stop_vars(1,1) = incrCost;                 obj.stop_vars(1,2) = 0;
            obj.stop_vars(2,1) = incrNormL2;              obj.stop_vars(2,2) = obj.max_constr_change;
            obj.stop_vars(3,1) = obj.line_search.kappa;     obj.stop_vars(3,2) = obj.line_search.kappa_min;
        end
        
    end
    
    methods (Access = public)
        
        function N_L2 = norm_L2(obj,x,x_ini)
            inc_x = x-x_ini;
            N_L2 = obj.scalar_product.computeSP_M(inc_x,inc_x)/obj.scalar_product.computeSP_M(x_ini,x_ini);
        end
        
    end
    
    methods
        
        function constr_tol = get.constr_tol(obj)
            constr_tol(1:obj.nconstr) = obj.target_parameters.constr_tol;
        end
        
    end
    
end