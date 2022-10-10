classdef DesignVariableAcademic < handle
    
    properties (Access = public)
        value
        valueOld
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = DesignVariableAcademic()
            
        end
        
        function init(obj,x0)
            obj.value = x0;
        end
        
    end
    
    methods (Access = public)
        
        function update(obj,x)
            obj.value = x;
        end
        
        function updateOld(obj)
            obj.valueOld = obj.value;
        end
       
    end
    
end