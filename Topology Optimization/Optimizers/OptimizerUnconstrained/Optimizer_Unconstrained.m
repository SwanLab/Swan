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
        xOld
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
            obj.objectiveFunction  = cParams.lagrangian;
            obj.hasConverged       = false;            
            obj.maxIncrNormX       = +Inf;            
            obj.line_search        = LineSearch.create(cParams.lineSearchSettings);
            obj.scalar_product     = ScalarProduct(cParams.scalarProductSettings);
            obj.convergenceVars    = cParams.convergenceVars;
            obj.target_parameters  = cParams.target_parameters;
            obj.designVariable     = cParams.designVariable;
        end
        
        function update(obj)
            obj.compute();
            obj.objectiveFunction.updateBecauseOfPrimal();
            obj.updateConvergenceParams();
            
            if ~obj.hasConverged
                obj.line_search.computeKappa();
            end
        end
        
        function init(obj)
            obj.objectiveFunction.setInitialValue();
            obj.initLineSearch();
            obj.hasConverged = false;            
        end
        

        
        function updateConvergenceParams(obj)
            incrementNormL2  = obj.norm_L2(obj.designVariable.value,obj.designVariable.valueOld);
            incrementObjFunc = obj.objectiveFunction.computeIncrement();
            
            obj.designImproved = incrementObjFunc < 0 && incrementNormL2 < obj.maxIncrNormX;
            
            obj.hasConverged = obj.designImproved || obj.line_search.kappa <= obj.line_search.kappa_min;
            
            obj.convergenceVars.reset();
            obj.convergenceVars.append(incrementObjFunc);
            obj.convergenceVars.append(incrementNormL2);
            obj.convergenceVars.append(obj.line_search.kappa);
        end
        
    end
    
    methods (Access = private)        
        
        function initLineSearch(obj)
            x0 = obj.designVariable.valueOld;
            g  = obj.objectiveFunction.gradient;
            obj.line_search.initKappa(x0,g);
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