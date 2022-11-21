classdef ISCSOMaterial < handle
    
    properties (Access = public)
        E
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ISCSOMaterial(cParams)
            obj.init(cParams)
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.E = 1;
        end
        
    end
    
end