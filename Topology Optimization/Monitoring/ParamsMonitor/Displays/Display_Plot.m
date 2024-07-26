classdef Display_Plot < Display_Abstract

    methods (Access = public)

        function obj = Display_Plot(cParams)
            obj@Display_Abstract(cParams.title)
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
                obj.valueArray(end+1,:) = value(1);
                if length(value) == 1
                    obj.iterationArray(end+1) = it;
                elseif length(value) == 2
                    obj.iterationArray(end+1) = value(2);
                else
                    error('Data input error')
                end
            end
        end

        function refresh(obj)
            if ~isempty(obj.valueArray) && ~isempty(obj.iterationArray)
                set(obj.handle,'XData',obj.iterationArray,'YData',obj.valueArray);
                if obj.iterationArray(end)>0
                    set(obj.style,'XLim',[min(0,min(obj.iterationArray)), max(obj.iterationArray)])
                end
                drawnow
            end
        end

    end
    
end