classdef TrussStructureCost < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunction
    end
    
    methods (Access = public)
        
        function obj = TrussStructureCost(cParams)
            f = ShapeFunctional_Factory();
            obj.shapeFunction = f.create(cParams);
        end

        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end

        function computeFunction(obj)
            obj.shapeFunction.computeFunction();
            obj.value = obj.shapeFunction.value;
        end

        function computeGradient(obj)
            obj.shapeFunction.computeGradient();
            obj.gradient = obj.shapeFunction.gradient;
        end
        
    end
    
    methods (Access = private)
        
        
    end
    
end