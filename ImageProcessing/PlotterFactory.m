classdef PlotterFactory < handle

    methods (Access = public, Static)

        function plotter = create(cParams)
            switch cParams.type
                case 'Density'
                    plotter   = PlotterDensity(cParams);
                case 'LevelSet'
                    plotter   = PlotterLevelSet(s);   
            end
        end

    end

end