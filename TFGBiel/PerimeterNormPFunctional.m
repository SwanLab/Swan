classdef PerimeterNormPFunctional < handle

    properties (Access = private)
        quadrature
        Pp
        totalVolume
        filter
        epsilon
    end

    properties (Access = private)
        mesh
        perimeterTarget
        p
    end

    methods (Access = public)
        function obj = PerimeterNormPFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
            obj.createFilter();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)   % dJ: integrand field
            [xD,Le] = obj.computeFilteredVariable(x); % rho
            J       = obj.computeFunction(xD{1},Le);
            dJ{1}   = obj.computeGradient(xD{1},Le);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh            = cParams.mesh;
            obj.perimeterTarget = cParams.perimeterTarget;
            obj.p               = cParams.p;
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
            obj.epsilon  = 6*obj.mesh.computeMeanCellSize();
            f.updateEpsilon(obj.epsilon);
            obj.filter = f;
        end

        function [xD,Le] = computeFilteredVariable(obj,x)
            xD = x.obtainDomainFunction();
            Le = obj.filter.compute(xD{1},3);
        end

        function J = computeFunction(obj,x,Le)
            xP     = (x.*(1-Le)).^obj.p;
            PerP   = Integrator.compute(xP,obj.mesh,obj.quadrature.order);
            obj.Pp = PerP^(1/obj.p);
            J      = ((1/obj.perimeterTarget)*((1/obj.totalVolume)^(1/obj.p))*obj.Pp) - 1;
        end

        function dJ = computeGradient(obj,x,Le)
            Lea = obj.computeFilteredTermForGradient(x,Le);
            num = (((x.*(1-Le)).^(obj.p-1)).*(1-Le) - Lea).*(obj.Pp^(1-obj.p));
            den = obj.perimeterTarget*(obj.totalVolume)^(1/obj.p);
            dJ  = num./den;
            dJ  = dJ.project('P1');
        end

        function Lea = computeFilteredTermForGradient(obj,x,Le)
            a = ((x.*(1-Le)).^(obj.p-1)).*x;
            Lea = obj.filter.compute(a,3);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Perimeter p-norm';
        end
    end
end