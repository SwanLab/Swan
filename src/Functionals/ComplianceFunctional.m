classdef ComplianceFunctional < handle

    properties (Access = private)
        value0
        oldCost
        oldGradient
        xOld
    end

    properties (Access = private)
        mesh
        filter
        compliance
        material
        plotTitle
    end

    methods (Access = public)
        function obj = ComplianceFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            dx = xR{1} - obj.xOld;
            if norm(dx.fValues)/norm(xR{1}.fValues) > 0.02
                obj.material.setDesignVariable(xR);
                [J,dJ] = obj.computeComplianceFunctionAndGradient(x);
                obj.oldCost = J;
                obj.oldGradient = dJ;
                obj.xOld = xR{1};
            else
                sp = ScalarProduct(obj.oldGradient{1},dx,'L2');
                J = obj.oldCost + sp;
                dJ = obj.oldGradient;
            end
        end

        function title = getTitleToPlot(obj)
            title = obj.plotTitle;
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.material   = cParams.material;
            obj.compliance = cParams.complainceFromConstitutive;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
            obj.xOld = 1000;

            if isfield(cParams, 'title')
                obj.plotTitle = cParams.title;
            else
                obj.plotTitle = 'Compliance';
            end
        end

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj,x)
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            dC  = ChainRule.compute(x,dC);
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
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

    % methods (Static, Access = public)
    %     function title = getTitleToPlot()
    %         title = 'Compliance';
    %     end
    % end
end