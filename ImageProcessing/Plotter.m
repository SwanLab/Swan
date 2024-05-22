classdef Plotter < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            obj = PlotterFactory.create(cParams);
        end

    end

end