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
        lineSearch
        scalar_product
        constr_tol
        convergenceVars
    end
    
    properties (Access = protected)                
        designVariable        
        xOld
    end
    
    properties (Access = private)
       incX
       incF
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
    
    methods (Access = protected)
        
        function opt = obtainOptimalityTolerance(obj)
            opt = obj.targetParameters.optimality_tol;            
        end

    end
    
    methods (Access = public)
        
        function obj = Optimizer_Unconstrained(cParams) 
            obj.objectiveFunction  = cParams.lagrangian;
            obj.hasConverged       = false;            
            obj.maxIncrNormX       = +Inf;            
            obj.lineSearch         = LineSearch.create(cParams.lineSearchSettings);
            obj.scalar_product     = ScalarProduct(cParams.scalarProductSettings);
            obj.convergenceVars    = cParams.convergenceVars;
            obj.targetParameters   = cParams.targetParameters;
            obj.designVariable     = cParams.designVariable;
        end
        
        function update(obj)
            obj.designVariable.updateOld();                                    
            obj.init();
            while ~obj.hasConverged
                obj.designVariable.restart();
                obj.compute();
                obj.objectiveFunction.updateBecauseOfPrimal();
                obj.updateConvergenceParams();
                if ~obj.hasConverged
                    obj.lineSearch.computeKappa();
               end
            end
            obj.revertIfDesignNotImproved();
        end
        
        function init(obj)
            obj.objectiveFunction.updateOld();
            obj.initLineSearch();
            obj.hasConverged = false;            
        end
        
        function storeKappaVariable(obj)
            obj.convergenceVars.set(obj.lineSearch.kappa);            
        end
        
        function updateConvergenceParams(obj)
            obj.computeIncrements();
            obj.computeOptimizerFlagConvergence(); 
            obj.storeConvergenceVariablesValues();            
        end
        
        function storeConvergenceVariablesValues(obj)
            obj.convergenceVars.reset();
            obj.convergenceVars.append(obj.incF);
            obj.convergenceVars.append(obj.incX);
            obj.convergenceVars.append(obj.lineSearch.kappa);            
        end
        
        function computeOptimizerFlagConvergence(obj)
            costDecreased = obj.hasCostDecreased();
            smallChangeX  = obj.isVariableChangeSmall();
            isValidIter   = costDecreased && smallChangeX;
            isLineSearchSmall = obj.isLineSearchSmallerThanMin();
            obj.hasConverged = isValidIter || isLineSearchSmall;            
        end
        
        function itHas = hasCostDecreased(obj)
            itHas = obj.incF < 0;
        end
        
        function itIs = isVariableChangeSmall(obj)
            itIs = obj.incX < obj.maxIncrNormX;
        end
        
        function itIs = isLineSearchSmallerThanMin(obj)
           itIs = obj.lineSearch.kappa <= obj.lineSearch.kappa_min;
        end
        
        function computeIncrements(obj)
            normXsquare = obj.designVariable.computeL2normIncrement();
            obj.incX = sqrt(normXsquare);
            obj.incF = obj.objectiveFunction.computeIncrement();
        end
        
        function itIs = isOptimal(obj)
            optimTol = obj.obtainOptimalityTolerance();
            optCond  = obj.opt_cond;
            isNot = optCond >= optimTol;
            itIs = ~isNot;
        end
        
       function initLineSearch(obj)
            x0 = obj.designVariable.value;
            g  = obj.objectiveFunction.gradient;
            obj.lineSearch.initKappa(x0,g);
        end
        
    end
    
    methods (Access = private)        
        

        
        function revertIfDesignNotImproved(obj)
            if ~obj.designImproved
                obj.designVariable.restart();
            end
        end
        
    end
    
%     methods
%         
%         function constr_tol = get.constr_tol(obj)
%             constr_tol = obj.targetParameters.constr_tol;
%         end
%         
%     end
%     
end