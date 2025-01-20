classdef VolumeNormPFunctional < handle

    properties (Access = private)
        quadrature
        totalVolume
        volume
    end

    properties (Access = private)
        mesh
        alpha
        p
        gradientTest
    end

    methods (Access = public)
        function obj = VolumeNormPFunctional(cParams)
            obj.init(cParams);
            obj.createGlobalVolume(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)   % dJ: integrand field
            xD    = x.obtainDomainFunction(); % rho
%             [V,~] = obj.volume.computeFunctionAndGradient(x);
            J     = obj.computeFunction(xD{1});
            dJ{1} = obj.computeGradient(xD,J);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.alpha = cParams.alpha;
            obj.p     = cParams.p;
            obj.gradientTest = cParams.gradientTest;
        end

        function createTotalVolume(obj)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            obj.totalVolume = sum(dV(:));
        end

        function createGlobalVolume(obj,cParams)
            obj.volume = VolumeFunctional(cParams);
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,3);
            obj.quadrature = quad;
        end

        function J = computeFunction(obj,x)
            xP   = x.^obj.p;
            volP = Integrator.compute(xP,obj.mesh,obj.quadrature.order);
            J    = volP^(1/obj.p);
            J    = J/obj.totalVolume - obj.alpha;
        end

        function dJ = computeGradient(obj,x,Vp)
            dx      = obj.gradientTest;
            rho     = x{1}.fValues;
            fValues = (Vp^(1-obj.p))*rho.^(obj.p-1);
            dJ      = FeFunction.create(dx.order,fValues,obj.mesh);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume p-norm';
        end
    end
end