classdef Display_SemilogY < Display_Abstract
    
    methods (Access = public)
        
        function obj = Display_SemilogY(title)
            obj@Display_Abstract(title);
        end
        
    end
   
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = semilogy(0,0);
        end
        
    end
    
end