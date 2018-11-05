classdef Tensor < handle
    
    properties (Access = protected)
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
    
    methods (Access = protected)
        
        function loadTensorVariables(obj)

        end
        
    end
    
    methods (Access = protected, Abstract)
%         loadOrderVariable(obj)
%         loadRepresentationVariable(obj)
%         loadElasticityCaseVariable(obj)        
    end
    
end

