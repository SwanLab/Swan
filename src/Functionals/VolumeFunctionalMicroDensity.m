classdef VolumeFunctionalMicroDensity < handle

    properties (Access = private)
        mesh
        base
        test
        baseFun
        totalVolume
    end

    methods (Access = public)

        function obj = VolumeFunctionalMicroDensity(cParams)
            obj.init(cParams);
            obj.createBaseFunction();
            obj.createTotalVolume();
        end

        function [J, dJ] = computeFunctionAndGradient(obj, x)
            xD    = x.obtainDomainFunction();
            rho   = xD{2};               % volume depende só de rho
            J     = obj.computeFunction(rho);
            dJ{1} = obj.computeGradientB();    % dV/db = 0
            dJ{2} = obj.computeGradientRho();  % dV/drho
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.base = cParams.uMesh;
            obj.test = cParams.test;
        end

        function createBaseFunction(obj)
            s.trial    = LagrangianFunction.create(obj.mesh, 1, obj.test.order);
            s.mesh     = obj.mesh;
            riszFilter = FilterLump(s);
            f          = CharacteristicFunction.create(obj.base);
            obj.baseFun = riszFilter.compute(f, 2);
        end

        function createTotalVolume(obj)
            V = Integrator.compute(obj.baseFun, obj.mesh, 2);
            obj.totalVolume = V;
        end

        function J = computeFunction(obj, rho)
            volume = Integrator.compute(rho .* obj.baseFun, obj.mesh, 2);
            J      = volume / obj.totalVolume;
        end

        function dJ = computeGradientRho(obj)
            % dV/drho = baseFun / totalVolume — igual ao original
            dJ = copy(obj.baseFun);
            dJ.setFValues(dJ.fValues ./ obj.totalVolume);
        end

        function dJ = computeGradientB(obj)
            % dV/db = 0 — volume nao depende de b
            dJ = copy(obj.baseFun);
            dJ.setFValues(zeros(size(obj.baseFun.fValues)));
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end

end
