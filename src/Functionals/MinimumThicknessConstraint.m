classdef MinimumThicknessConstraint < handle

    properties (Access = private)
        functional
        hmin
        hmin0
        valueOld
        tarVolume
        totalVolume
    end

    methods (Access = public)
        function obj = MinimumThicknessConstraint(cParams)
            obj.init(cParams);
            obj.computeTotalVolume(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            V  = Integrator.compute(xD{1},x.fun.mesh,2);
            V  = V/(obj.tarVolume*obj.totalVolume) - 1;
            obj.updateTarget(x.fun,V);

            [Jf,dJf] = obj.functional.computeFunctionAndGradient(x);
            J = 1-Jf/obj.hmin0;
            dJ{1} = (-1/obj.hmin0).*dJf{1};

            J = max(J,-1); % To avoid bad scaling of monitoring
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.functional = MinimumThicknessFunctional(cParams);
            obj.hmin = cParams.target;
            obj.hmin0 = cParams.target0;
            obj.valueOld = -inf;
            obj.tarVolume = cParams.tarVolume;
        end

        function computeTotalVolume(obj,s)
            u = ConstantFunction.create(1,s.mesh);
            obj.totalVolume = Integrator.compute(u,s.mesh,2);
        end

        function updateTarget(obj,x,V) % Cuando la suma de grays empieza a decaer puede provocar tmb la decay
            if norm(x.fValues-obj.valueOld)==0 && V<0.1
                obj.hmin0 = min(obj.hmin0/0.8,obj.hmin);
            end
            obj.valueOld = x.fValues;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Min thick constr';
        end
    end
end