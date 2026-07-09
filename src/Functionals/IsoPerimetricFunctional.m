classdef IsoPerimetricFunctional < handle

    properties (Access = private)
        mesh
        base
        epsilon
        p
    end

    properties (Access = private)
        filter
        riszFilter
        baseFun
        totalVolume
    end

    methods (Access = public)
        function obj = IsoPerimetricFunctional(cParams)
            obj.init(cParams);
            obj.createFilter();
            obj.createRiszFilter();
            obj.createBaseFunction();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [xD,Le] = obj.computeFilteredVariable(x);
            Qp      = obj.computeIsoPerimetricQuocient(xD{1},Le);
            J       = obj.computeFunction(Qp);
            dJ{1}   = obj.computeGradient(xD{1},Le,Qp);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.base    = cParams.uMesh;
            obj.epsilon = cParams.epsilon;
            obj.p       = cParams.p;
        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            f.updateEpsilon(obj.epsilon);
            obj.filter = f;
        end

        function createRiszFilter(obj)
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            obj.riszFilter = FilterLump(s);
        end

        function createBaseFunction(obj)
            f           = CharacteristicFunction.create(obj.base);
            obj.baseFun = obj.riszFilter.compute(f,2);
        end

        function createTotalVolume(obj)
            dV = obj.baseFun;
            V  = Integrator.compute(dV,obj.mesh,2);
            obj.totalVolume = V;
        end

        function [xD,Le] = computeFilteredVariable(obj,x)
            xD = x.obtainDomainFunction();
            Le = obj.filter.compute(xD{1},3);
        end

        function Qp = computeIsoPerimetricQuocient(obj,x,Le)
            b       = obj.baseFun;
            xP      = (((x.^1.5).*(1-Le).*((x.*Le).^(-0.5))).*(1/(2*obj.epsilon))).^obj.p;
            isoPerP = Integrator.compute(xP.*b,obj.mesh,3);
            Qp      = isoPerP^(1/obj.p);
        end

        function J = computeFunction(obj,Qp)
            J = (1/obj.totalVolume)^(1/obj.p)*Qp;
        end

        function dJ = computeGradient(obj,x,Le,Qp)
            bChi = obj.baseFun;
            b    = ((x.^1.5).*(1-Le).*((x.*Le).^(-0.5))).^(obj.p-1);
            Lea  = computeFilteredTermForGradient(obj,x,Le,b);
            num1 = (b.^1.5).*((b.*Le).^(-0.5) - (b^(-1)).*(b.*Le).^(0.5)).*(Qp.^(1-obj.p));
            num2 = Lea.*(Qp.^(1-obj.p));
            den  = ((2*obj.epsilon)^obj.p)*(obj.totalVolume)^(1/obj.p);
            dJ1  = num1./den;
            dJ1  = obj.riszFilter.compute(dJ1.*bChi,3);
            dJ2  = num2./den;
            dJ2  = obj.riszFilter.compute(dJ2.*bChi,3);
            dJ   = dJ1 + dJ2;
        end

        function Lea = computeFilteredTermForGradient(obj,x,Le,b)
            a   = (x.^2.5).*(((x.*Le).^(-1.5))./(-2) - (x^(-1)).*(x.*Le).^(-0.5)./(-2)).*b;
            Lea = obj.filter.compute(a,3);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'IsoPerimeter fun';
        end
    end
end