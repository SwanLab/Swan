classdef BoundayConditionsPlotter_Null < BoundayConditionsPlotter_Abstract
    
    methods (Access = public)
        
        function obj = BoundayConditionsPlotter_Null()
        end
        
        function plotDirichlet(~)
        end
        
        function plotNeumann(~)
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj)
        end
        
    end
    
end