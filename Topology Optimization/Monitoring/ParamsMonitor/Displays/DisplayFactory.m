classdef DisplayFactory < handle
    
    methods (Access = public, Static)
        
        function display = create(s)
            switch s.chartType
                case 'plot'
                    display = Display_Plot(s);
                case 'log'
                    display = Display_SemilogY(s);
                case 'bar'
                    display = Display_Bar(s);
                case 'surf'
                    display = Display_Surf(s);
                otherwise
                    error('Invalid Chart Type.')
            end
        end
        
    end
    
end