classdef FlexureDOFConstraint < handle

    properties (Access = private)
        strainEnergyFunc
        gscale
    end

    methods (Access = public)
        function obj = FlexureDOFConstraint(cParams)
            obj.strainEnergyFunc = cParams.strainEnergyFunc;
            obj.gscale = cParams.gscale;
        end

        function [J, dJ] = computeFunctionAndGradient(obj,x)
            [E,dE] = obj.strainEnergyFunc.computeFunctionAndGradient(x);
            J = obj.gscale*(E-1);

            for i = 1:length(dE)
                dE{i}.setFValues(dE{i}.fValues*obj.gscale);
            end
            dJ = dE;
        end

        function title = getTitleToPlot(obj)
            title = 'DOF Flexibility Constraint';
        end
    end
end

