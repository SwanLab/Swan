classdef CostComplianceVolume < handle

    properties (Access = private)
        lambda
        compliance
        volume
    end

    properties (Access = private)
        designVariable
        topOpt
    end

    methods (Access = public)

        function obj = CostComplianceVolume(cParams)
            obj.init(cParams);
            obj.createCost();
        end

        function [J,dJ] = computeValueAndGradient(obj,x)
            obj.designVariable.update(x);
            obj.compliance.computeFunctionAndGradient();
            obj.volume.computeFunctionAndGradient();
            J  = obj.computeCost;
            dJ = obj.computeGradient;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.topOpt         = cParams.topOpt;
            obj.designVariable = cParams.designVariable;
            obj.lambda = 5;
        end

        function J = computeCost(obj)
            c = obj.compliance.value;
            v = obj.volume.value;
            l = obj.lambda;
            J = c + l*v;
        end

        function dJ = computeGradient(obj)
            dc = obj.compliance.gradient;
            dv = obj.volume.gradient;
            l = obj.lambda;
            dJ = dc + l*dv;
        end

        function createCost(obj)
            obj.compliance = obj.topOpt.cost.shapeFunctions{1};
            obj.volume = obj.topOpt.constraint.shapeFunctions{1};5;
        end

    end

end