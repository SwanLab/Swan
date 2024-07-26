classdef Display_Abstract < handle
    
    properties (Access = protected)
        handle
        iterationArray
        valueArray
        style

        figTitle
    end

    methods (Access = protected, Abstract)
        setChartType(obj)
    end

    methods (Access = public , Static)
        
        function obj = create(s)
            d = DisplayFactory();
            obj = d.create(s);
        end
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
        
        function updateParams(obj,it,value)
            obj.iterationArray(end+1) = it;
            if ~isempty(value)
                obj.valueArray(end+1,:) = value;
            end
        end

        function refresh(obj)
            if ~isempty(obj.valueArray) && ~isempty(obj.iterationArray)
                set(obj.handle,'XData',obj.iterationArray,'YData',obj.valueArray);
                if obj.iterationArray(end)>0
                    set(obj.style,'XLim',[0 obj.iterationArray(end)])
                end
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

