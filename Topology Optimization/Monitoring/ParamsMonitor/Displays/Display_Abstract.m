classdef Display_Abstract < handle
    
    properties (Access = protected)
        handle
        iterationArray
        valueArray
        style
    end
    
    properties (Access = private)
        figTitle
    end
    
    methods (Access = protected, Abstract)
        
        setChartType(obj)
        
    end
    
    methods (Access = public)
        
        function obj = Display_Abstract(figTitle)
            obj.figTitle = figTitle;
        end
        
        function show(obj,nRows,nCols,i,margins)
            obj.setStyle(nRows,nCols,i,margins);
            obj.setChartType();
            title(obj.figTitle);
            grid on
        end
        
        function updateParams(obj,it,value,monitor)
            obj.iterationArray(end+1) = it;
            if ~isempty(value)
                obj.valueArray(end+1,:) = value;
            end
            monitor.increaseRefreshIterator();
        end
        
        function refresh(obj)
            if ~isempty(obj.valueArray) && ~isempty(obj.iterationArray)
                set(obj.handle,'XData',obj.iterationArray,'YData',obj.valueArray);
                set(obj.style,'XLim',[0 obj.iterationArray(end)])
                drawnow
            end
        end
        
    end
    
    methods (Access = private)
        
        function setStyle(obj,nRows,nCols,i,margins)
            obj.style = subplot_tight(nRows,nCols,i,margins);
        end
        
    end
    
end

