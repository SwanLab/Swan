classdef Display_Plot < Display_Abstract
    
    methods (Access = public)
        
        function obj = Display_Plot(title)
            obj@Display_Abstract(title);
        end
        
    end
   
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = plot(0,0);
        end
        
    end
    
end