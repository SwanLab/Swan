classdef DisplaySemilogX < DisplayPlot
    
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = semilogx(0,0);
        end
        
    end
    
end