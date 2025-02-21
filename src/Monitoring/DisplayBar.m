classdef DisplayBar < DisplayPlot
    
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = bar(0,0);
        end
        
    end

    methods (Access = public)

        function refresh(obj)
            if ~isempty(obj.ArrayDataX) && ~isempty(obj.ArrayDataY)
                axes  = obj.obtainDisplayAxes();
                nTick = 5;
                nP    = length(obj.ArrayDataX);
                if nP<=nTick
                    set(axes,'XTick',obj.ArrayDataX);
                else
                    step = ceil(nP/nTick);
                    set(axes,'XTick',obj.ArrayDataX(1:step:end));
                end
                set(obj.handle,'XData',obj.ArrayDataX,'YData',obj.ArrayDataY);
                if obj.ArrayDataY(end)>0
                    set(obj.style,'XLim',[min(0,min(obj.ArrayDataX)), max(0,max(obj.ArrayDataX)+1)])
                    set(obj.style,'YLim',[min(-0.01,min(obj.ArrayDataY)), max(0,max(obj.ArrayDataY)+1)])
                end
            end
            grid off;
        end

    end
    
end