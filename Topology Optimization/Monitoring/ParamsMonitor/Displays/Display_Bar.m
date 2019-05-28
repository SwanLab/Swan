classdef Display_Bar < Display_Abstract
    
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = bar(0,0);
        end
        
    end
    
end