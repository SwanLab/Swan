classdef VolumeFunctional < handle

    properties (Access = private)
        quadrature
        totalVolume
    end

    properties (Access = private)
        mesh
        gradientTest
    end

    methods (Access = public)
        function obj = VolumeFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD    = x.obtainDomainFunction();
            J     = obj.computeFunction(xD);
            dJ{1} = obj.computeGradient(xD);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.gradientTest = cParams.gradientTest;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createTotalVolume(obj)
            vt = obj.mesh.computeVolume();
            obj.totalVolume = vt;
        end

        function J = computeFunction(obj,x)
            volume = Integrator.compute(x,obj.mesh,obj.quadrature.order);
            J      = volume/obj.totalVolume;
        end

        function dJ = computeGradient(obj,x)
            dJ = obj.gradientTest.copy();
            dJ.fValues = ones(dJ.nDofs,1)/obj.totalVolume;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end

