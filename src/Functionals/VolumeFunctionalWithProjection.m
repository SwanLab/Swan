classdef VolumeFunctionalWithProjection < handle

    properties (Access = private)
        quadrature
        totalVolume
    end

    properties (Access = private)
        mesh
        gradientTest
        iter
        filter
        filterAdjoint
        filteredDesignVariable
    end

    methods (Access = public)
        function obj = VolumeFunctionalWithProjection(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
%             iter = x{2};
%             x = x{1};

            xD  = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            J  = obj.computeFunction(xR{1});
            dJ = obj.computeGradient(xR);
            obj.filteredDesignVariable = xR{1};
        end

        function x = getDesignVariable(obj)
            x = obj.filteredDesignVariable;
        end
        
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.gradientTest = cParams.gradientTest;
            obj.filter       = cParams.filter;
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint = cParams.filterAdjoint;            
            end
            obj.iter = 0;
        end

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
                if ~isempty(obj.filterAdjoint)
                    obj.filterAdjoint.updateFilteredField(obj.filter);
                end
            end
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end

        function createTotalVolume(obj)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            obj.totalVolume = sum(dV(:));
        end

        function J = computeFunction(obj,x)
            volume = Integrator.compute(x,obj.mesh,obj.quadrature.order);
            J      = volume/obj.totalVolume;
        end

        function dJ = computeGradient(obj,x)
            test    = obj.gradientTest;
            fValues = ones(test.nDofs,1)/obj.totalVolume;
            dJ      = FeFunction.create(test.order,fValues,obj.mesh);
            if ~isempty(obj.filterAdjoint)
                dJ     = obj.filterAdjoint.compute(dJ,2);
            else
                dJ     = obj.filter.compute(dJ,2);
            end
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end

