classdef IsotropicPerimeterNormPFunctional < handle

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
        function obj = IsotropicPerimeterNormPFunctional(cParams)
            obj.init(cParams);
            obj.createFilter();
            obj.createRiszFilter();
            obj.createBaseFunction();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [xD,Le] = obj.computeFilteredVariable(x);
            Pp      = obj.computePerimeterPNorm(xD{1},Le);
            J       = obj.computeFunction(Pp);
            dJ{1}   = obj.computeGradient(xD{1},Le,Pp);
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

        function Pp = computePerimeterPNorm(obj,x,Le)
            b    = obj.baseFun;
            xP   = ((x.*(1-Le)).*(1/(2*obj.epsilon))).^obj.p;
            PerP = Integrator.compute(xP.*b,obj.mesh,3);
            Pp   = PerP^(1/obj.p);
        end

        function J = computeFunction(obj,Pp)
            J = ((1/obj.totalVolume)^(1/obj.p))*Pp;
        end

        function dJ = computeGradient(obj,x,Le,Pp)
            b   = obj.baseFun;
            Lea = obj.computeFilteredTermForGradient(x,Le);
            num = ((((x.*(1-Le)).*(1/(2*obj.epsilon))).^(obj.p-1)).*(1-Le) - Lea).*(Pp^(1-obj.p));
            den = 2*obj.epsilon*(obj.totalVolume)^(1/obj.p);
            dJ  = num./den;
            dJ  = obj.riszFilter.compute(dJ.*b,3);
        end

        function Lea = computeFilteredTermForGradient(obj,x,Le)
            a = x.*(((x.*(1-Le)).*(1/(2*obj.epsilon))).^(obj.p-1));
            Lea = obj.filter.compute(a,3);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Iso Per p-norm';
        end
    end
end