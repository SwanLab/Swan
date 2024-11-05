classdef DisplayAbstract < handle
    
    properties (Access = protected)
        figTitle

        position
        handle
        iteration
        ArrayDataX
        ArrayDataY
        style
    end

    methods (Access = protected, Abstract)
        setChartType(obj)
    end

    methods (Access = public, Abstract)
        updateParams(ibj,it,value)
    end

    methods (Access = public , Static)
        
        function obj = create(s)
            d = DisplayFactory();
            obj = d.create(s);
        end
    end

    methods (Access = public)

        function obj = DisplayAbstract(figTitle,position)
            obj.figTitle = figTitle;
            obj.position = position;
        end

        function show(obj,nRows,nCols,i,margins)
            obj.setStyle(nRows,nCols,i,margins);
            obj.setChartType();
            title(obj.figTitle);
            grid on
        end

        function refresh(obj)
            if ~isempty(obj.valueArray) && ~isempty(obj.iterationArray)
                set(obj.handle,'XData',obj.ArrayDataX,'YData',obj.ArrayDataY);
                if obj.ArrayDataX(end)>0
                    set(obj.style,'XLim',[0 obj.ArrayDataX(end)])
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


