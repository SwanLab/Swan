classdef Optimizer_Unconstrained < handle
    
    properties (Access = public)
        objectiveFunction
        targetParameters
        opt_cond
        hasConverged        
    end
    
    properties (GetAccess = public, SetAccess = protected)
        designImproved
        maxIncrNormX
    end
    
    properties (GetAccess = public, SetAccess = private)
        line_search
        scalar_product
        constr_tol
        convergenceVars
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
            obj.targetParameters   = cParams.targetParameters;
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
            nIncX = obj.designVariable.computeL2normIncrement();
            nIncF = obj.objectiveFunction.computeIncrement();
            
            obj.designImproved = nIncF < 0 && nIncX < obj.maxIncrNormX;
            isLineSearchSmallerThanMin = obj.line_search.kappa <= obj.line_search.kappa_min;
            
            obj.hasConverged = obj.designImproved || isLineSearchSmallerThanMin;
            
            obj.convergenceVars.reset();
            obj.convergenceVars.append(nIncF);
            obj.convergenceVars.append(nIncX);
            obj.convergenceVars.append(obj.line_search.kappa);
        end
        
    end
    
    methods (Access = private)        
        
        function initLineSearch(obj)
            x0 = obj.designVariable.value;
            g  = obj.objectiveFunction.gradient;
            obj.line_search.initKappa(x0,g);
        end
        
    end
    
    methods
        
        function constr_tol = get.constr_tol(obj)
            constr_tol = obj.targetParameters.constr_tol;
        end
        
    end
    
end