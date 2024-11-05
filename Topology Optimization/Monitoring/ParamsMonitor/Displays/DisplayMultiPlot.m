classdef DisplayMultiPlot < DisplayAbstract
    
    properties (Access = public)
        legend
    end
    methods (Access = public)

        function obj = DisplayMultiPlot(cParams)
            obj@DisplayAbstract(cParams.title,cParams.position)
            obj.legend = cParams.legend;
        end
        
    end


    methods (Access = protected)
        
        function setChartType(obj)
            nLines = length(obj.legend);
            for i=1:nLines
                obj.handle{i} = plot(1,1);
                hold on
            end
            legend(obj.legend);
        end
        
    end

    methods (Access = public)

        function updateParams(obj,it,value)

            if ~isempty(value)
                nLines = length(obj.legend);
                obj.ArrayDataY(end+1,:) = value(1:nLines);
                if length(value) == nLines
                    obj.ArrayDataX(end+1) = it;
                elseif length(value) == nLines+1
                    obj.ArrayDataX(end+1) = value(end);
                else
                    error('Data input error')
                end
            end
        end

        function refresh(obj)
            if ~isempty(obj.ArrayDataX) && ~isempty(obj.ArrayDataY)
                nLines = length(obj.legend);
                for i=1:nLines
                    set(obj.handle{i},'XData',obj.ArrayDataX,'YData',obj.ArrayDataY(:,i));
                    if obj.ArrayDataX(end)>0
                        set(obj.style,'XLim',[min(0,min(obj.ArrayDataX)), max(obj.ArrayDataX)])    
                    end
                    drawnow
                end
            end
        end

    end
    
end