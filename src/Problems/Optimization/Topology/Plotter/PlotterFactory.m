classdef PlotterFactory < handle

    methods (Access = public, Static)

        function plotter = create(designVariable)
            d = designVariable;

            switch d.type
                case {'Density','DensityAndBound'}
                    s.mesh           = d.fun.mesh;
                    s.designVariable = d;
                    plotter          = PlotterDensity(s);

                case 'LevelSet'
                    s.designVariable = d;
                    plotter          = PlotterLevelSet(s);

                case 'MultiField'
                    s.designVariable = d;
                    plotter          = PlotterMultiField(s);

                otherwise
                    error('PlotterFactory:UnknownType', ...
                        'Unknown design variable type: %s',d.type);
            end
        end

    end

end