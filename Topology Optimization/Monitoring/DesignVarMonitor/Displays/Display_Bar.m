classdef Display_Bar < Display_Abstract
    
    methods (Access = public)
        
        function obj = Display_Bar(title)
            obj@Display_Abstract(title);
        end
        
    end
    
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = bar(0,0);
        end
        
    end
    
end