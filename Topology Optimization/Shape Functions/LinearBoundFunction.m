classdef LinearBoundFunction < handle

    methods (Static, Access = public)
        function [J,dJ] = computeFunctionAndGradient(x)
            J               = x.bound;
            dJ.fValues      = zeros(length(x.fun.fValues),1);
            dJ.fValues(end) = 1;
            dJ              = {dJ};
        end

        function title = getTitleToPlot()
            title = 'LinearBoundFun';
        end
    end

end