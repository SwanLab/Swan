classdef ComplianceFunctional < handle
    properties (Access = private)
        value0
    end
    properties (Access = private)
        mesh
        filter
        compliance
        C
        dC
    end
    methods (Access = public)
        function obj = ComplianceFunctional(cParams)
            obj.init(cParams);
        end
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            Crho = obj.C(xR{1});
            dCrho = obj.dC(xR{1});
            [J,dJ] = obj.computeComplianceFunctionAndGradient(x,Crho,dCrho);
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.C          = cParams.C;
            obj.dC         = cParams.dC;
            obj.compliance = cParams.complainceFromConstitutive;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end
        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end
        function [J,dJ] = computeComplianceFunctionAndGradient(obj,x,C,dRhoC)
            dxC{1} = ChainRule.compute(x,dRhoC);
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dxC);
            dJ     = obj.filterFields(dJ);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J  = obj.computeNonDimensionalValue(J);
            dJ = obj.computeNonDimensionalGradient(dJ);
        end
        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end
        function dx = computeNonDimensionalGradient(obj,dx)
            refX = obj.value0;
            for i = 1:length(dx)
                dx{i}.setFValues(dx{i}.fValues/refX);
            end
        end
    end
    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end