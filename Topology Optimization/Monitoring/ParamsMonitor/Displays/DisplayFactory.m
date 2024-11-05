classdef DisplayFactory < handle
    
    methods (Access = public, Static)
        
        function display = create(s)
            switch s.chartType
                case 'plot'
                    display = DisplayPlot(s);
                case 'multiPlot'
                    display = DisplayMultiPlot(s);
                case 'log'
                    display = Display_SemilogY(s);
                case 'bar'
                    display = Display_Bar(s);
                case 'surf'
                    display = DisplaySurf(s);
                otherwise
                    error('Invalid Chart Type.')
            end
        end
        
    end
    
end