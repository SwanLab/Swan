classdef DisplayLogLog < DisplayPlot
    
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = loglog(0,0);
        end
        
    end
    
end