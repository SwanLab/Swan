classdef Tensor < handle
    
    properties (Abstract)
        tensor
    end
    
    methods (Access = public)
        
        function T = getValue(obj)
           T = obj.tensor; 
        end
        
        function setValue(obj,T)
            obj.tensor = T;
        end
    end
    
end

