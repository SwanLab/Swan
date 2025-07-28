classdef VolumeFunctional < handle

    properties (Access = private)
        mesh
        base
        test
    end

    properties (Access = private)
        riszFilter
        baseFun
        totalVolume
    end

    methods (Access = public)
        function obj = VolumeFunctional(cParams)
            obj.init(cParams);
            obj.createRiszFilter();
            obj.createBaseFunction();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            J  = obj.computeFunction(xD{1});
            dJ{1} = obj.computeGradient();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.base = cParams.uMesh;
            obj.test = cParams.test;
        end

        function createRiszFilter(obj)
            s.trial = LagrangianFunction.create(obj.mesh,1,obj.test.order);
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

        function J = computeFunction(obj,x)
            b      = obj.baseFun;
            volume = Integrator.compute(x.*b,obj.mesh,2);
            J      = volume/obj.totalVolume;
        end

        function dJ = computeGradient(obj)
            dj = obj.baseFun./obj.totalVolume;
            dJ = obj.riszFilter.compute(dj,2);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume'; % Maybe a property in the future?
        end
    end
end

