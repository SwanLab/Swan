classdef Objective_Function < handle
    
    properties
        value
        gradient
    end
    
    properties (Access = protected)
        valueInitial
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
