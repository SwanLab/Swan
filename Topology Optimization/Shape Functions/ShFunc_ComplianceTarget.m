classdef ShFunc_ComplianceTarget < handle

    properties (Access = public)
        value
        gradient
        Msmooth
        dvolu
        value0
        shFunTargettedType
        shapeFunctionTargetted
    end

    properties (Access = protected)
        designVariable
    end

    methods (Access = public)

        function obj = ShFunc_ComplianceTarget(cParams)
            obj.createBound(cParams);
            obj.createFunctionalWithTarget(cParams);
            obj.setUpMassMatrixAndVolumen();
            obj.setUpInitialValue();
        end

        function computeFunctionAndGradient(obj)
            obj.shapeFunctionTargetted.computeFunctionAndGradient();
            obj.value = obj.shapeFunctionTargetted.value - obj.designVariable.bound;
            obj.gradient = [obj.shapeFunctionTargetted.gradient;-1];
            if isempty(obj.value0)
                obj.value0 = obj.shapeFunctionTargetted.value0;
            end 
        end

        function f = getPhysicalProblems(obj)
            f = obj.shapeFunctionTargetted.getPhysicalProblems();
        end

        function f = getRegularizedDesignVariable(obj)
            f = obj.shapeFunctionTargetted.getRegularizedDesignVariable();
        end

        function q = getQuad(obj)
            q = obj.shapeFunctionTargetted.getQuad();
        end
        function fP = addPrintableVariables(obj)
            fP = obj.shapeFunctionTargetted.addPrintableVariables();
        end
        function v = getVariablesToPlot(obj)
            v = obj.shapeFunctionTargetted.getVariablesToPlot();
        end
        function t = getTitlesToPlot(obj)
            t = obj.shapeFunctionTargetted.getTitlesToPlot();
        end
        function fP = createPrintVariables(obj)
            fP = obj.shapeFunctionTargetted.createPrintVariables();
        end
        function [fun, funNames] = getFunsToPlot(obj)
             [fun, funNames] = obj.shapeFunctionTargetted.getFunsToPlot();
        end
        function updateTargetParameters(obj)
            obj.shapeFunctionTargetted.updateTargetParameters();
        end

    end
    methods (Access = protected)
        
        function setUpMassMatrixAndVolumen(obj)
            obj.Msmooth = obj.shapeFunctionTargetted.Msmooth;
            obj.dvolu = obj.shapeFunctionTargetted.dvolu;
        end 
        function setUpInitialValue(obj)
            obj.value0 = obj.shapeFunctionTargetted.value0;
        end 
        function createFunctionalWithTarget(obj,cParams)
            obj.shFunTargettedType = cParams.filterParams.femSettings.shFunType;
            cParams.type  = obj.shFunTargettedType;
            cParams.designVariable = cParams.designVariable.density;
            f = ShapeFunctional_Factory();
            obj.shapeFunctionTargetted = f.create(cParams);
        end
        function createBound(obj,cParams)
            obj.designVariable = cParams.designVariable;
        end
    end
end