classdef Damage < DesignVariable

    methods (Access = public)

        function obj = Damage(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
        end

        function fxV = evaluate(obj,xV)
            fxV = obj.fun.evaluate(xV);
        end

        function grad = computeGrad(obj)
            grad = obj.fun.computeGrad();
        end

        function update(obj,value)
            obj.fun.setFValues(value);
        end

        function plot(obj)
            plot(obj.fun);
        end
    
    end
    
end

