classdef GradientMethod < handle
    
    properties (Access = private)
        lineSearch
        differentiableFunction
        designVariable
    end

    methods (Access = public)
        
        function obj = GradientMethod(cParams)
            obj.differentiableFunction = cParams.differentiableFunction;
            obj.designVariable = cParams.designVariable;
            obj.createLineSearch()
        end

        function compute(obj)
            t = obj.lineSearch.value;
            obj.differentiableFunction.computeGradient();
            g = obj.differentiableFunction.gradient;
            x = obj.designVariable.value;
            x = x - t*g;
            obj.designVariable.value = x;
        end

    end
    
    methods (Access = private)
        
        function createLineSearch(obj)
            L = obj.differentiableFunction.lipschitzConstant;
            s.value = 1/L;
            obj.lineSearch = ConstantLineSearch(s);
        end
        
    end
    
end