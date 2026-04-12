classdef MinimumThicknessFunctional < handle

    properties (Access = private)
        tau
        volume
        perimeter
        baseFun
        totalVolume
    end

    methods (Access = public)
        function obj = MinimumThicknessFunctional(cParams)
            obj.init(cParams);
            obj.createBaseFunction(cParams);
            obj.computeTotalVolume(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            Vo = obj.totalVolume;
            [V,dV] = obj.volume.computeFunctionAndGradient(x);
            V = V*Vo;
            dV{1}.setFValues(dV{1}.fValues.*Vo);

            [P,dP] = obj.perimeter.computeFunctionAndGradient(x);

            t = obj.tau;
            f = Vo*t^2/(V+t^2);

            J = 2*(V+f)/(P+t);
            dJ{1} = (2/(P+t)).*(1-Vo*t^2/(V+t^2)^2).*dV{1} - (2*(V+f)/((P+t)^2)).*dP{1};

            s.trial     = LagrangianFunction.create(x.fun.mesh,1,'P1');
            s.mesh = x.fun.mesh;
            riszFilter  = FilterLump(s);
            dJ{1} = riszFilter.compute(dJ{1},2);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.tau       = cParams.tau;
            obj.volume    = VolumeFunctional(cParams);
            obj.perimeter = PerimeterFunctional(cParams);
        end

        function createBaseFunction(obj,s)
            s.trial     = LagrangianFunction.create(s.mesh,1,'P1');
            riszFilter  = FilterLump(s);
            f           = CharacteristicFunction.create(s.uMesh);
            obj.baseFun = riszFilter.compute(f,2);
        end

        function computeTotalVolume(obj,s)
            dV = obj.baseFun;
            V  = Integrator.compute(dV,s.mesh,2);
            obj.totalVolume = V;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Thickness';
        end
    end
end