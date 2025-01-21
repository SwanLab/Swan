classdef VolumeNormPFunctional < handle

    properties (Access = private)
        quadrature
        totalVolume
        filter
    end

    properties (Access = private)
        mesh
        alpha
        p
    end

    methods (Access = public)
        function obj = VolumeNormPFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
            obj.createFilter();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)   % dJ: integrand field
            xD    = x.obtainDomainFunction(); % rho
            J     = obj.computeFunction(xD{1});
            dJ{1} = obj.computeGradient(xD,J);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.alpha = cParams.alpha;
            obj.p     = cParams.p;
        end

        function createTotalVolume(obj)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            obj.totalVolume = sum(dV(:));
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,3);
            obj.quadrature = quad;
        end

        function J = computeFunction(obj,x)
            xP   = x.^obj.p;
            volP = Integrator.compute(xP,obj.mesh,obj.quadrature.order);
            J    = (volP/obj.totalVolume)^(1/obj.p);
            J    = J/obj.alpha - 1;
        end

        function dJ = computeGradient(obj,x,Vp)
            rho = x{1};
            dj  = ((Vp^(1-obj.p))*rho.^(obj.p-1))./(obj.alpha*obj.totalVolume^(1/obj.p));
            dJ  = obj.filter.compute(dj,3);
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume p-norm';
        end
    end
end