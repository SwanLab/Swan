classdef NumericalGradientComputer < handle
    properties (Access = public)
        gradient
    end
    properties (Access = private)
        xC
        s
        objective
        objBase
        xp
        bp
    end

    methods (Access = public)
        function obj = NumericalGradientComputer(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.computeObjectiveBase();
            obj.computeGradient();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.xC = cParams.x;
            obj.s = cParams.s;
            obj.bp = cParams.bp;
        end

        function computeObjectiveBase(obj)
            u.bp = obj.bp;
            u.x = obj.xC;
            objB = ObjectiveFunctionComputer(u);
            objB.compute();
            obj.objBase = objB.objectiveFunc;
        end

        function computeObjective(obj)
            u.bp = obj.bp;
            u.x = obj.xp;
            object = ObjectiveFunctionComputer(u);
            object.compute();
            obj.objective = object.objectiveFunc;
        end

        function computeGradient(obj)
            x = obj.xC;
            n = size(x,2);
            m = size(obj.s,2);
            ep = 1e-5;
            for i = 1:n
                obj.xp = x;
                obj.xp(i) = x(i) + ep;
                obj.computeObjective();
                g(i) = [(obj.objective - obj.objBase)/ep]';
            end
            if (m >= 1)
                g(n+1:n+m) = 0;
            end
            obj.gradient = g;
        end
    end
end