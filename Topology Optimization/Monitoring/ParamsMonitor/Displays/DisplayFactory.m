classdef DisplayFactory < handle
    
    methods (Access = public, Static)
        
        function display = create(s)
            switch s.chartType
                case 'plot'
                    display = DisplayPlot(s);
                case 'multiplot'
                    display = DisplayMultiPlot(s);
                case 'semilogY'
                    display = DisplaySemilogY(s);
                case 'loglog'
                    display = DisplayLogLog(s);
                case 'bar'
                    display = DisplayBar(s);
                case 'surf'
                    display = DisplaySurf(s);
                otherwise
                    error('Invalid Chart Type.')
            end
        end
        
    end
    
end