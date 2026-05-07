classdef FlexureDOFConstraint < handle

    properties (Access = private)
        strainEnergyFunc
        gscale
        emax
    end

    methods (Access = public)
        function obj = FlexureDOFConstraint(cParams)
            obj.strainEnergyFunc = cParams.strainEnergyFunc;
            obj.gscale = cParams.gscale;
            obj.emax = cParams.emax;
        end

        function [J, dJ] = computeFunctionAndGradient(obj,x)
            [E,dE] = obj.strainEnergyFunc.computeFunctionAndGradient(x);
            value0 = obj.strainEnergyFunc.value0;
            J = obj.gscale*(E*value0-obj.emax);
            %J = obj.gscale*(E-1);

            for i = 1:length(dE)
                dE{i}.setFValues(dE{i}.fValues*value0*obj.gscale);
            end
            dJ = dE;
        end

        function title = getTitleToPlot(obj)
            title = 'DOF Flexibility Constraint';
        end
    end
end

