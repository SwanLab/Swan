classdef IsoPerimetricNormPFunctional < handle

    properties (Access = private)
        quadrature
        Qp
        totalVolume
        filterDesignVariable
        filterAdjoint
        epsilon
    end

    properties (Access = private)
        mesh
        C
        p
    end

    methods (Access = public)
        function obj = IsoPerimetricNormPFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
            obj.createFilter();
            obj.createFilterAdjoint();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [xD,Le] = obj.computeFilteredVariable(x);
            J       = obj.computeFunction(xD{1},Le);
            dJ{1}   = obj.computeGradient(xD{1},Le);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.C    = cParams.C;
            obj.p    = cParams.p;
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
            obj.filterDesignVariable = f;
        end

        function createFilterAdjoint(obj)
            s.filterType = 'PDE';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterAdjoint = f;
        end

        function [xD,Le] = computeFilteredVariable(obj,x)
            xD = x.obtainDomainFunction();
            Le = obj.filterDesignVariable.compute(xD{1},3);
        end

        function J = computeFunction(obj,x,Le)
            xP      = ((x.*(1-Le).*(Le.^(-0.5))).*(1/(2*obj.epsilon))).^obj.p;
            IsoPerP = Integrator.compute(xP,obj.mesh,obj.quadrature.order);
            obj.Qp  = IsoPerP^(1/obj.p);
            J       = ((1/obj.C)*((1/obj.totalVolume)^(1/obj.p))*obj.Qp) - 1;
        end

        function dJ = computeGradient(obj,x,Le)
            b   = (x.*(1-Le).*(Le.^(-0.5))).^(obj.p-1);
            Lea = computeFilteredTermForGradient(obj,x,Le,b);
            num = (obj.Qp.^(1-obj.p)).*(((Le.^(-0.5) - Le.^(0.5)).*b) + Lea);
            den = ((2*obj.epsilon)^obj.p)*obj.C*(obj.totalVolume)^(1/obj.p);
            dJ  = num./den;
            dJ  = obj.filterAdjoint.compute(dJ,3);
        end

        function Lea = computeFilteredTermForGradient(obj,x,Le,b)
            a   = x.*(-0.5.*Le.^(-1.5) - 0.5.*Le.^(-0.5)).*b;
            Lea = obj.filterDesignVariable.compute(a,3);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'IsoPerimeter p-norm';
        end
    end
end