classdef ObjectiveFunction < handle
    
    properties (Access = public)
        value
        gradient
    end
    
    properties (Access = protected)
        cost
        constraint  
        dualVariable
    end
    
    properties (Access = private)
        valueOld        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.cost         = cParams.cost;
            obj.constraint   = cParams.constraint;   
            obj.dualVariable = cParams.dualVariable;
        end
        
    end
    
    methods (Access = public)
        
        function incr = computeIncrement(obj)
            v  = obj.value;
            vI = obj.valueOld;
            incr = (v - vI)/abs(vI);
        end
        
        function setInitialValue(obj)
            obj.valueOld = obj.value;
        end
        
        function restart(obj)
            obj.value = obj.valueOld;
        end
        
        
    end
    
end
