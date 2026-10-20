classdef ComplianceFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filters
        compliance
        C
        dC
    end

    methods (Access = public)

        function obj = ComplianceFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            Cx = obj.C(xR);
            dCx = obj.dC(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient(x,Cx,dCx);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.C          = cParams.C;
            obj.dC         = cParams.dC;
            obj.compliance = cParams.complainceFromConstitutive;

            if isfield(cParams,'filters')
                obj.filters = cParams.filters;
            elseif isfield(cParams,'filter')
                obj.filters = {cParams.filter};
            else
                error('ComplianceFunctional:MissingFilter', ...
                    'A filter or a cell array of filters must be provided.');
            end

            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end

        function xR = filterFields(obj,x)
            nDesVar = numel(x);
            nFilters = numel(obj.filters);

            if nFilters ~= 1 && nFilters ~= nDesVar
                error('ComplianceFunctional:InvalidNumberOfFilters', ...
                    'Expected either one filter or one filter per design variable.');
            end

            xR = cell(size(x));

            for i = 1:nDesVar
                if nFilters == 1
                    filterI = obj.filters{1};
                else
                    filterI = obj.filters{i};
                end

                if isempty(filterI)
                    xR{i} = x{i};
                else
                    xR{i} = filterI.compute(x{i},2);
                end
            end
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj,x,Cx,dCx)
            dxC = ChainRule.compute(x,dCx);
            [J,dJ] = obj.compliance.computeFunctionAndGradient(Cx,dxC);
            dJ = obj.filterFields(dJ);

            if isempty(obj.value0)
                obj.value0 = J;
            end

            J = obj.computeNonDimensionalValue(J);
            dJ = obj.computeNonDimensionalGradient(dJ);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x = x/refX;
        end

        function dx = computeNonDimensionalGradient(obj,dx)
            refX = obj.value0;

            for i = 1:numel(dx)
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