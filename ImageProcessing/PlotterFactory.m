classdef PlotterFactory < handle

    methods (Access = public, Static)

        function plotter = create(cParams)
            switch cParams.type
                case 'Density'
                    s.density = cParams.density;
                    plotter   = PlotterDensity(s);
                case 'LevelSet'
                    s.designVariable = d;
                    plotter          = PlotterLevelSet(s);   
            end
        end

    end

end