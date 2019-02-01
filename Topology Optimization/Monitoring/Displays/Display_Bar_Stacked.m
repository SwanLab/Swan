classdef Display_Bar_Stacked < Display_Bar
    
    methods (Access = public)
        
        function obj = Display_Bar_Stacked(title)
            obj@Display_Bar(title);
        end
        
        function refresh(obj)
            if ~isempty(obj.valueArray) && ~isempty(obj.iterationArray)                
                it = obj.iterationArray(end);
                value = obj.valueArray(end,:);
                set(obj.handle,'XData',0,'YData',value(1));
                set(obj.style,'XLim',[-1 1],'YLim',[0 value(2)])
                set(gca, 'XTick', []);
                drawnow
            end
        end
        
    end
    
end