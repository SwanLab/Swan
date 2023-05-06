classdef AcademicCost < handle
    
    properties (Access = private)
        costFunction
        gradientFunction
        designVariable
    end
    
    properties (Access = public)
        value
        gradient
    end
    
    methods (Access = public)
        
        function obj = AcademicCost(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj)
            x            = obj.designVariable.value';
            obj.value    = obj.costFunction(x);
            obj.gradient = obj.gradientFunction(x);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.costFunction     = cParams.cF;
            obj.gradientFunction = cParams.gF;
            obj.designVariable   = cParams.dV;
        end
        
    end
    
end