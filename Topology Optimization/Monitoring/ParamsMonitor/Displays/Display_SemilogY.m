classdef Display_SemilogY < Display_Plot
    
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = semilogy(0,0);
        end
        
    end
    
end