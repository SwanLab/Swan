classdef VolumeFunctional < handle

    properties (Access = private)
        quadrature
        totalVolume
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        function obj = VolumeFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            J  = obj.computeFunction(xD);
            dJ = obj.computeGradient(xD);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createTotalVolume(obj)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            obj.totalVolume = sum(dV(:));
        end

        function J = computeFunction(obj,x)
            int    = Integrator.create('Function',obj.mesh,obj.quadrature.order);
            volume = int.compute(x);
            J      = volume/obj.totalVolume;
        end

        function dJ = computeGradient(obj,x)
            fValues = ones(x.nDofs,1)/obj.totalVolume;
            dJ      = FeFunction.create(x.order,fValues,obj.mesh);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end

