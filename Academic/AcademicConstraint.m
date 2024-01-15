classdef AcademicConstraint < handle
    
    properties (Access = private)
        gradientFunction
        constraintFunction
        designVariable
    end
    
    properties (Access = public)
        value
        gradient
        nSF
    end
    
    methods (Access = public)
        
        function obj = AcademicConstraint(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj)
            x            = obj.designVariable.value;
            obj.value    = obj.constraintFunction(x);
            obj.gradient = obj.gradientFunction(x);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.constraintFunction = cParams.cF;
            obj.gradientFunction   = cParams.gF;
            obj.designVariable     = cParams.dV;
            obj.nSF                = cParams.nSF;
        end
        
    end
    
end