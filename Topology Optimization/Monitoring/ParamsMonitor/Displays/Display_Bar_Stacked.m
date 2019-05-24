classdef Display_Bar_Stacked < Display_Bar
    
    methods (Access = public)
        
        function refresh(obj)
            if ~isempty(obj.valueArray) && ~isempty(obj.iterationArray)                
                it = obj.iterationArray(end);
                value = obj.valueArray(end,:);
                set(obj.handle,'XData',0,'YData',value);
                set(obj.style,'XLim',[-1 1],'YLim',[0 1])
                set(gca, 'XTick', []);
                drawnow
            end
        end
        
    end
    
end