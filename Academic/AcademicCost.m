classdef AcademicCost < handle
    
    properties (Access = private)
        costFunction
        gradientFunction
        designVariable
        hessianFunction
    end
    
    properties (Access = public)
        value
        gradient
        hessian
    end
    
    methods (Access = public)
        
        function obj = AcademicCost(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj)
            x            = obj.designVariable.value';
            obj.value    = obj.costFunction(x);
            obj.gradient = obj.gradientFunction(x);
            %obj.hessian = obj.hessianFunction(x);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.costFunction     = cParams.cF;
            obj.gradientFunction = cParams.gF;
            obj.designVariable   = cParams.dV;
            %obj.hessianFunction = cParams.hF;
        end
        
    end
    
end