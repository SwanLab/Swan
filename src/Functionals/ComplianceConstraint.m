classdef ComplianceConstraint < handle

    properties (Access = private)
        mesh
        complianceTarget
        compliance
    end

    methods (Access = public)

        function obj = ComplianceConstraint(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [c,dc] = obj.compliance.computeFunctionAndGradient(x);
            J      = obj.computeFunction(c);
            dJ     = obj.computeGradient(dc);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.complianceTarget = cParams.complianceTarget;
            obj.compliance       = ComplianceFunctional(cParams);
        end

        function J = computeFunction(obj,c)
            cTar = obj.complianceTarget;
            J    = c/cTar-1;
        end

        function dJ = computeGradient(obj,dc)
            cTar    = obj.complianceTarget;
            nDesVar = numel(dc);
            dJ      = cell(size(dc));

            for i = 1:nDesVar
                fValues = dc{i}.fValues/cTar;
                dJ{i}   = FeFunction.create(dc{i}.order,fValues,obj.mesh);
            end
        end

    end

    methods (Static, Access = public)

        function title = getTitleToPlot()
            title = 'Compliance constraint';
        end

    end

end