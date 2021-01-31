classdef DisplayFactory < handle
    
    methods (Access = public, Static)
        
        function display = create(chartType,title)
            switch chartType
                case 'plot'
                    display = Display_Plot(title);
                case 'log'
                    display = Display_SemilogY(title);
                case 'bar'
                    display = Display_Bar(title);           
                otherwise
                    error('Invalid Chart Type.')
            end
        end
        
    end
    
end