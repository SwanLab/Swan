classdef VolumeFunctionalMicroDensity < handle

    properties (Access = private)
        mesh
        base
        test
    end

    properties (Access = private)
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
            xD = x.obtainDomainFunction();
            
            % CORREÇÃO: usar a segunda variável (ρ) para o volume
            rho = xD{2};  % ← mudar de {1} para {2}
            
            J = obj.computeFunction(rho);
            
            % Gradiente para ρ (segunda variável)
            dJ_rho = obj.computeGradient();
            
            % Gradiente para b (primeira variável) é ZERO
            dJ_b = dJ_rho.copy();
            dJ_b.setFValues(zeros(size(dJ_b.fValues)));
            
            % Retornar ambos os gradientes
            dJ = {dJ_b, dJ_rho};
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.base = cParams.uMesh;
            obj.test = cParams.test;
        end
        
        function createBaseFunction(obj)
            s.trial     = LagrangianFunction.create(obj.mesh, 1, obj.test.order);
            s.mesh      = obj.mesh;
            riszFilter  = FilterLump(s);
            f           = CharacteristicFunction.create(obj.base);
            obj.baseFun = riszFilter.compute(f, 2);
        end

        function createTotalVolume(obj)
            dV = obj.baseFun;
            V = Integrator.compute(dV, obj.mesh, 2);
            obj.totalVolume = V;
        end

        function J = computeFunction(obj, rho)
            b = obj.baseFun;
            volume = Integrator.compute(rho .* b, obj.mesh, 2);
            J = volume / obj.totalVolume;
        end

        function dJ = computeGradient(obj)
            dJ = copy(obj.baseFun);
            dJ.setFValues(dJ.fValues ./ obj.totalVolume);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end