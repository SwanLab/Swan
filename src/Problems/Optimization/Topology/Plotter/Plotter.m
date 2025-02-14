classdef Plotter < handle

    methods (Static, Access = public)

        function obj = create(type)
            obj = PlotterFactory.create(type);
        end

    end

end