classdef AreaColumn < DesignVariable
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = AreaColumn(cParams)
            obj.init(cParams)
            
        end
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value;
        end        
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            
        end
        
    end
    
end