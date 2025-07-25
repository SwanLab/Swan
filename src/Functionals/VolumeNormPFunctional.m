classdef VolumeNormPFunctional < handle

    properties (Access = private)
        quadrature
        vP
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
            xD    = obj.computeFilteredVariable(x); % rho
            J     = obj.computeFunction(xD);
            dJ{1} = obj.computeGradient(xD);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.alpha = cParams.alpha;
            obj.p     = cParams.p;
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,3);
            obj.quadrature = quad;
        end

        function createTotalVolume(obj)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            obj.totalVolume = sum(dV(:));
        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            epsilon      = 6*obj.mesh.computeMeanCellSize();
            f.updateEpsilon(epsilon);
            obj.filter = f;
        end

        function xD = computeFilteredVariable(obj,x)
            xD = x.obtainDomainFunction();
            xD = obj.filter.compute(xD{1},3);
        end

        function J = computeFunction(obj,x)
            xP     = x.^obj.p;
            volP   = Integrator.compute(xP,obj.mesh,obj.quadrature.order);
            obj.vP = volP^(1/obj.p);
            J      = (volP/obj.totalVolume)^(1/obj.p);
            J      = J/obj.alpha - 1;
        end

        function dJ = computeGradient(obj,x)
            rho = x;
            dj  = ((obj.vP^(1-obj.p))*rho.^(obj.p-1))./(obj.alpha*obj.totalVolume^(1/obj.p));
            dJ  = obj.filter.compute(dj,3);
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume p-norm';
        end
    end
end