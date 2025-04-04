classdef ShFunc_Compliance_constraint < ShFunWithElasticPdes
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        shapeFunctionCompliance
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Compliance_constraint(cParams)
            obj.init(cParams);
            obj.constraintValue = cParams.constraintValue;
            obj.createShapeFunctionCompliance(cParams);
        end
        
        function computeFunctionAndGradient(obj)
            obj.shapeFunctionCompliance.computeFunctionAndGradient();
        end
        
        function computeFunction(obj)
            obj.shapeFunctionCompliance.computeFunction();
            obj.value = obj.shapeFunctionCompliance.value - obj.constraintValue;
        end
        
        function computeGradient(obj)
            obj.shapeFunctionCompliance.computeGradient();
            obj.gradient = obj.shapeFunctionCompliance.gradient;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
        function createShapeFunctionCompliance(cParams)
            obj.shapeFunctionCompliance = ShFunc_Compliance(cParams);
        end
        
    end
    
end