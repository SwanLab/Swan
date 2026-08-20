classdef MotionBasedStiffnessConstraint < handle

    % Constraint on the motion based characteristic stiffness of a free MP.
    % g = gscale ( Kp(x)/Kp_bar - 1)<=0
    % with Kp = 2*E/gamma^2, which for unit displ in gamma is Kp=2*E

    properties (Access = private)
        MotionBasedStrainEnergy
        gscale
        Kp_bar
    end

    methods (Access = public)
        function obj = MotionBasedStiffnessConstraint(cParams)
            obj.MotionBasedStrainEnergy = cParams.MotionBasedStrainEnergy;
            obj.gscale = cParams.gscale;
            obj.Kp_bar = cParams.Kp_bar;
        end

        function [J, dJ] = computeFunctionAndGradient(obj,x)
            [E,dE] = obj.MotionBasedStrainEnergy.computeFunctionAndGradient(x);
            value0 = obj.MotionBasedStrainEnergy.value0;
            J = obj.gscale*(2*E*value0/obj.Kp_bar - 1);

            for i = 1:length(dE)
                dE{i}.setFValues(dE{i}.fValues*value0*obj.gscale*2/obj.Kp_bar);
            end
            dJ = dE;
        end

        function title = getTitleToPlot(obj)
            title = 'Motion Based Stiffness Constraint';
        end
    end
end

