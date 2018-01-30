classdef ShFunc_ChomogLamPerimeter_alphabeta< Shape_Functional
    properties
        chomog
        perimeter
        lambda
    end
    methods
        function obj=ShFunc_ChomogLamPerimeter_alphabeta(settings)
            obj@Shape_Functional(settings);
            obj.perimeter=ShFunc_Perimeter(settings);
            obj.chomog=ShFunc_Chomog_alphabeta(settings);
            obj.lambda=settings.perimeter.lambda;
        end
    end
    methods
        function computef(obj,x,physicalProblem,interpolation,filter)            
            obj.perimeter.computef(x, physicalProblem, interpolation,filter);
            obj.chomog.computef(x, physicalProblem, interpolation,filter);
            
            obj.value=obj.chomog.value+obj.lambda*obj.perimeter.value;
            obj.gradient=obj.chomog.gradient+obj.lambda*obj.perimeter.gradient;
        end
    end
end