classdef DesignVariable < handle
    
    properties (Access = public)
       value
       mesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = DesignVariable()
            obj.init()            
        end
        
        function update(obj,x)
            obj.value = x;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
    end
    
end