classdef MicroDamageParams < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = MicroDamageParams(cParams)
            obj.nVariables = 1;
            obj.init(cParams)
            
        end
        
        function update(obj,x)
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end