classdef VolumeConstraintWithBound < handle

    properties (Access = private)
        volume
    end

    methods (Access = public)
        function obj = VolumeConstraintWithBound(cParams)
            obj.volume = VolumeConstraint(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD         = x.density;
            [J,dJv]    = obj.volume.computeFunctionAndGradient(xD);
            dJ.fValues = [dJv{1}.fValues;0];
            dJ         = {dJ};
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume constraint';
        end
    end
end