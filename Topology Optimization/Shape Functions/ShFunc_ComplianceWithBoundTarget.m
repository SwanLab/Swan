classdef ShFunc_ComplianceWithBoundTarget < handle

    properties (Access = public)
        value
        gradient
        shFunTargettedType
    end

    properties (Access = private)
        shapeFunctionTargetted
        designVariable
    end

    methods (Access = public)

        function obj = ShFunc_ComplianceWithBoundTarget(cParams)
            obj.init(cParams);
            obj.createFunctionalWithTarget(cParams);
        end

        function computeFunctionAndGradient(obj)
            obj.shapeFunctionTargetted.computeFunctionAndGradient();
            obj.value = obj.shapeFunctionTargetted.value - obj.designVariable.bound;
            obj.gradient = [obj.shapeFunctionTargetted.gradient;-1];
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

        function updateTargetParameters(obj)
            obj.shapeFunctionTargetted.updateTargetParameters();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
        end

        function createFunctionalWithTarget(obj,cParams)
            obj.shFunTargettedType = cParams.filterParams.femSettings.shFunType;
            cParams.type  = obj.shFunTargettedType;
            cParams.designVariable = cParams.designVariable.density;
            f = ShapeFunctional_Factory();
            obj.shapeFunctionTargetted = f.create(cParams);
        end

    end
end