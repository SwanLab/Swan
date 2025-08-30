classdef VolumeFunctional < handle

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
        function obj = VolumeFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
%             iter = x{2};
%             x = x{1};
            xD  = x.obtainDomainFunction();
            J  = obj.computeFunction(xD{1});
            dJ = obj.computeGradient(xD);

            obj.filteredDesignVariable = xD{1};
        end

        function x = getDesignVariable(obj)
            x = obj.filteredDesignVariable;
        end
        
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.gradientTest = cParams.gradientTest;
            obj.iter = 0;
            if isfield(cParams,'filterAdjoint')
                obj.filter        = cParams.filter;
                obj.filterAdjoint = cParams.filterAdjoint;
            end
        end

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
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
            end
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end

