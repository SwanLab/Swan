classdef BoundayConditionsPlotter_Null < BoundayConditionsPlotter_Abstract
    
    methods (Access = public)
        
        function obj = BoundayConditionsPlotter_Null(axes,mesh)
            obj@BoundayConditionsPlotter_Abstract(axes,mesh);
        end
        
        function plotDirichlet(obj)
            
        end
        
        function plotNeumann(obj)
            
        end
        
    end
    
end