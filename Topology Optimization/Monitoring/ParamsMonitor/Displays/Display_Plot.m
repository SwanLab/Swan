classdef Display_Plot < Display_Abstract

    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = plot(0,0);
        end
        
    end
    
end