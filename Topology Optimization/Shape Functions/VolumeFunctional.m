classdef VolumeFunctional < handle

    properties (Access = private)
        quadrature
        totalVolume
    end

    properties (Access = private)
        mesh
        gradientTest
        filter
        filterAdjoint
        filteredDesignVariable
        iter
    end

    methods (Access = public)
        function obj = VolumeFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            iter = x{2};
            x = x{1};
 
%             if iter > obj.iter
%                 obj.iter = iter;
%                 beta = obj.filter.getBeta();
%                 if iter >= 400 && mod(iter,20)== 0 && beta <= 40
%                     obj.filter.updateBeta(beta+2.0);
%                     obj.filterAdjoint.updateBeta(beta+2.0);
%                 end
%             end    

            xD  = x.obtainDomainFunction();
            xR  = obj.filterDesignVariable(xD);
            obj.filteredDesignVariable = xR;
            J  = obj.computeFunction(xR);
            dJ = obj.computeGradient(xR);
        end

        function dV = getDesignVariable(obj)
            dV = obj.filteredDesignVariable;
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

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
            if ~isempty(obj.filterAdjoint)
                xFiltered = obj.filter.onlyFilter(x,2);
                obj.filterAdjoint.updateFilteredField(xFiltered);
            end
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end

