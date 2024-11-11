classdef DisplaySemilogY < DisplayPlot
    
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = semilogy(0,0);
        end
        
    end
    
end