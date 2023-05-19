classdef ShFunTargetComputer < handle
    properties (Access = public)
        value
        gradient
        filter
        Msmooth
        dvolu
        value0
        shNumber
        regDesignVariable
        target
    end
    properties (Access = protected)
        ShapeFunction
    end
    methods (Access = public)
        function obj = ShFunTargetComputer(cParams)
            obj.createShFunc(cParams);
            obj.setUpMassMatrixAndVolumen();
            obj.setUpInitialValue();
        end
        function computeFunctionAndGradient(obj)
            obj.ShapeFunction.computeFunctionAndGradient();
            obj.value = obj.ShapeFunction.value/obj.target - 1;
            obj.gradient = obj.ShapeFunction.gradient;
        end
        function computeFunction(obj)
            obj.ShapeFunction.computeFunction();
            obj.value = obj.ShapeFunction.value/obj.target- 1;
        end

        function computeGradient(obj)
            obj.ShapeFunction.computeGradient();
            obj.gradient = obj.ShapeFunction.gradient;
        end

        function f = getPhysicalProblems(obj)
            f = obj.ShapeFunction.getPhysicalProblems(),
        end

        function f = getRegularizedDesignVariable(obj)
            f = obj.ShapeFunction.getRegularizedDesignVariable();
        end

        function q = getQuad(obj)
            q = obj.ShapeFunction.getQuad();
        end
        function fP = addPrintableVariables(obj)
            fP = obj.ShapeFunction.addPrintableVariables();
        end
        function v = getVariablesToPlot(obj)
            v = obj.ShapeFunction.getVariablesToPlot();
        end
        function t = getTitlesToPlot(obj)
            t = obj.ShapeFunction.getTitlesToPlot();
        end
        function fP = createPrintVariables(obj)
            fP = obj.ShapeFunction.createPrintVariables();
        end
        function [fun, funNames] = getFunsToPlot(obj)
             [fun, funNames] = obj.ShapeFunction.getFunsToPlot();
        end
        function updateTargetParameters(obj)
            obj.ShapeFunction.updateTargetParameters();
        end

    end
    methods (Access = protected)
        function setUpMassMatrixAndVolumen(obj)
            obj.Msmooth = obj.ShapeFunction.Msmooth;
            obj.dvolu = obj.ShapeFunction.dvolu;
        end 
        function setUpInitialValue(obj)
            obj.value0 = obj.ShapeFunction.value0;
        end 
    end
end