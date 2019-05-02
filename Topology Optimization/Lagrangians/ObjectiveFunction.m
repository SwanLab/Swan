classdef ObjectiveFunction < handle
    
    properties (Access = public)
        value
        gradient
    end
    
    properties (Access = protected)
        valueInitial
        cost
        constraint  
        dualVariable
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
            vI = obj.valueInitial;
            incr = (v - vI)/abs(vI);
        end
        
        function setInitialValue(obj)
            obj.valueInitial = obj.value;
        end
        
    end
    
end
