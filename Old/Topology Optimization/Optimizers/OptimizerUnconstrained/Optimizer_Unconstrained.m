classdef Optimizer_Unconstrained < handle

    properties (Access = public)
        objectiveFunction
        targetParameters
        optimalityCond
        hasConverged
    end

    properties (GetAccess = public, SetAccess = protected)
        designImproved
        maxIncrNormX
        minIncrNormX
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
        incX
        incF
    end

    properties (Access = private)

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



    end

    methods (Access = public)

        function obj = Optimizer_Unconstrained(cParams)
            obj.objectiveFunction  = cParams.lagrangian;
            obj.hasConverged       = false;
            obj.minIncrNormX       = 1e-8;%1e-2;
            obj.maxIncrNormX       = 5*1e-2;
            obj.convergenceVars    = cParams.convergenceVars;
            obj.targetParameters   = cParams.targetParameters;
            obj.designVariable     = cParams.designVariable;
            obj.createScalarProductCalculator(cParams);
            obj.createLineSearch(cParams);
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
                    obj.updateLineSearch();
                end
            end
            obj.revertIfDesignNotImproved();
        end

        function init(obj)
            obj.objectiveFunction.updateOld();
            obj.tryLineSearch();
            obj.hasConverged = false;
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
            obj.convergenceVars.append(obj.lineSearch.value);
            obj.convergenceVars.append(obj.lineSearch.nTrials);
        end

        function computeOptimizerFlagConvergence(obj)
            costDecreased    = obj.hasCostDecreased();
            xNotLargyChanged = obj.hasXNotLarglyChanged();
            isIterValid   = costDecreased && xNotLargyChanged;
            isLSsmall = obj.isLineSearchTooSmall();
            obj.hasConverged = isIterValid || isLSsmall;
        end

        function itHas = hasCostDecreased(obj)
            itHas = obj.incF < 0;
        end

        function itIs = isVariableChangeSmall(obj)
            itIs = obj.incX < obj.minIncrNormX;
        end

        function itIs = hasXNotLarglyChanged(obj)
            itIs = obj.incX < obj.maxIncrNormX;
        end
        
        function computeIncrements(obj)
            normXsquare = obj.designVariable.computeL2normIncrement();
            obj.incX = sqrt(normXsquare);
            obj.incF = obj.objectiveFunction.computeIncrement();
        end

        function itIs = isOptimal(obj)
            optimTol = obj.obtainOptimalityTolerance();
            optCond  = obj.optimalityCond;
            isNot = optCond >= optimTol;
            smallChangeX  = obj.isVariableChangeSmall();
            itIs = ~isNot || smallChangeX ;
        end

        function tryLineSearch(obj)
            obj.lineSearch.computeTrial();
        end
        
        function startLineSearch(obj)
            obj.lineSearch.computeStartingValue();
        end
        
        function itIs = isLineSearchTooSmall(obj)
            itIs = obj.lineSearch.isTooSmall();
        end
        
        function updateLineSearch(obj)
            obj.lineSearch.update();
        end

    end

    methods (Access = private)
        
        function opt = obtainOptimalityTolerance(obj)
            opt = obj.targetParameters.optimality_tol;
        end
        
        function createScalarProductCalculator(obj,cParams)
            s = cParams.scalarProductSettings;
            if isa(s,'ScalarProduct')
                obj.scalar_product = s;
            else
                obj.scalar_product = ScalarProduct(s);
            end
        end
        
        function createLineSearch(obj,cParams)
            s = cParams.lineSearchSettings;
            s.lineSearchInitiatorSettings.scalarProduct = obj.scalar_product;
            s.designVariable    = obj.designVariable;
            s.objectiveFunction = obj.objectiveFunction;
            obj.lineSearch  = LineSearch.create(s);
        end

        function revertIfDesignNotImproved(obj)
            if ~obj.designImproved
                obj.designVariable.restart();
            end
        end

    end

end
