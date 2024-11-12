classdef DisplayPlot < DisplayAbstract

    methods (Access = public)

        function obj = DisplayPlot(cParams)
            obj@DisplayAbstract(cParams)
        end
        
    end

    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = plot(0,0);
        end
        
    end

    methods (Access = public)

        function updateParams(obj,it,value)
            if ~isempty(value)
                obj.ArrayDataY(end+1,:) = value(1);
                if length(value) == 1
                    obj.ArrayDataX(end+1) = it;
                elseif length(value) == 2
                    obj.ArrayDataX(end+1) = value(2);
                else
                    error('Data input error')
                end
            end
        end

        function refresh(obj)
            if ~isempty(obj.ArrayDataX) && ~isempty(obj.ArrayDataY)
                set(obj.handle,'XData',obj.ArrayDataX,'YData',obj.ArrayDataY);
                if obj.ArrayDataY(end)>0
                    set(obj.style,'XLim',[min(0,min(obj.ArrayDataX)), max(1e-15,max(obj.ArrayDataX))])
                end
                drawnow;
            end
        end

    end
    
end