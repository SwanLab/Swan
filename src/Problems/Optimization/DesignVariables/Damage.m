classdef Damage < DesignVariable

    properties (Access = public)
        ndimf
        mesh
    end

    methods (Access = public)

        function obj = Damage(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
            obj.ndimf = obj.fun.ndimf;
            obj.mesh = obj.fun.mesh;
        end

        function fxV = evaluate(obj,xV)
            fxV = obj.fun.evaluate(xV);
        end

        function grad = computeGrad(obj,xV)
            grad = obj.fun.computeGrad(xV);
        end

        function update(obj,value)
            obj.fun.setFValues(value);
        end

        function plot(obj)
            plot(obj.fun);
        end
    
    end
    
end

