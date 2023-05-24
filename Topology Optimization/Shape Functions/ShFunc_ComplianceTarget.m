classdef ShFunc_ComplianceTarget < handle

    properties (Access = public)
        value
        gradient
        filter
        Msmooth
        dvolu
        value0
        shNumber
        regDesignVariable
    end

    properties (Access = protected)
        complianceShFunc
        designVariable
    end

    methods (Access = public)

        function obj = ShFunc_ComplianceTarget(cParams)
            obj.createBound(cParams);
            obj.createCompliance(cParams);
            obj.setUpMassMatrixAndVolumen();
            obj.setUpInitialValue();
        end

        function computeFunctionAndGradient(obj)
            obj.complianceShFunc.computeFunctionAndGradient();
            obj.value = obj.complianceShFunc.value - obj.designVariable.bound;
            obj.gradient = [obj.complianceShFunc.gradient;-1];
            if isempty(obj.value0)
                obj.value0 = obj.complianceShFunc.value0;
            end 
        end

        function f = getPhysicalProblems(obj)
            f = obj.complianceShFunc.getPhysicalProblems();
        end

        function f = getRegularizedDesignVariable(obj)
            f = obj.complianceShFunc.getRegularizedDesignVariable();
        end

        function q = getQuad(obj)
            q = obj.complianceShFunc.getQuad();
        end
        function fP = addPrintableVariables(obj)
            fP = obj.complianceShFunc.addPrintableVariables();
        end
        function v = getVariablesToPlot(obj)
            v = obj.complianceShFunc.getVariablesToPlot();
        end
        function t = getTitlesToPlot(obj)
            t = obj.complianceShFunc.getTitlesToPlot();
        end
        function fP = createPrintVariables(obj)
            fP = obj.complianceShFunc.createPrintVariables();
        end
        function [fun, funNames] = getFunsToPlot(obj)
             [fun, funNames] = obj.complianceShFunc.getFunsToPlot();
        end
        function updateTargetParameters(obj)
            obj.complianceShFunc.updateTargetParameters();
        end

    end
    methods (Access = protected)
        
        function setUpMassMatrixAndVolumen(obj)
            obj.Msmooth = obj.complianceShFunc.Msmooth;
            obj.dvolu = obj.complianceShFunc.dvolu;
        end 
        function setUpInitialValue(obj)
            obj.value0 = obj.complianceShFunc.value0;
        end 
        function createCompliance(obj,cParams)
            cParams.type  = 'compliance';
            cParams.designVariable = cParams.designVariable.density;
            f = ShapeFunctional_Factory();
            obj.complianceShFunc = f.create(cParams);
        end
        function createBound(obj,cParams)
            obj.designVariable = cParams.designVariable;
        end
    end
end