classdef Optimizer_Unconstrained < Optimizer
    
    properties (Access = public)
        objectiveFunction
        opt_cond
    end
    
    properties (GetAccess = public, SetAccess = protected)
        designImproved
        maxIncrNormX
    end
    
    properties (GetAccess = public, SetAccess = private)
        line_search
        scalar_product
        constr_tol
    end
    
    properties (Access = protected)
        designVariable        
    end
    
    methods (Access = public, Abstract)
        compute(obj)
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = Optimizer_UnconstrainedFactory();
            obj = f.create(cParams);
        end
        
        
    end
    
    methods (Access = public)
        
        function obj = Optimizer_Unconstrained(cParams)            
            obj.hasConverged   = false;            
            obj.maxIncrNormX   = +Inf;            
            obj.line_search    = LineSearch.create(cParams.lineSearchSettings);
            obj.scalar_product = ScalarProduct(cParams.scalarProductSettings);
            obj.convergenceVars = ConvergenceVariables(3);
            obj.target_parameters = cParams.target_parameters;
            obj.designVariable = cParams.designVariable;
        end
        
        function x = update(obj,x0)
            x = obj.compute(x0,obj.objectiveFunction.gradient);
            obj.designVariable.value = x;
            obj.objectiveFunction.updateBecauseOfPrimal(x);
            obj.updateConvergenceParams(x,x0);
            
            if ~obj.hasConverged
                obj.line_search.computeKappa();
            end
        end
        
        function init(obj,x0,objFunc)
            obj.objectiveFunction = objFunc;
            obj.objectiveFunction.setInitialValue();
            obj.initLineSearch(x0);
            obj.hasConverged = false;
        end
        
    end
    
    methods (Access = private)
        
        function updateConvergenceParams(obj,x,x_ini)
            incrementNormL2  = obj.norm_L2(x,x_ini);
            incrementObjFunc = obj.objectiveFunction.computeIncrement();
            
            obj.designImproved = incrementObjFunc < 0 && incrementNormL2 < obj.maxIncrNormX;
            
            obj.hasConverged = obj.designImproved || obj.line_search.kappa <= obj.line_search.kappa_min;
            
            obj.convergenceVars.reset();
            obj.convergenceVars.append(incrementObjFunc);
            obj.convergenceVars.append(incrementNormL2);
            obj.convergenceVars.append(obj.line_search.kappa);
        end
        
        function initLineSearch(obj,x0)
            obj.line_search.initKappa(x0,obj.objectiveFunction.gradient);
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
            constr_tol = obj.target_parameters.constr_tol;
        end
        
    end
    
end