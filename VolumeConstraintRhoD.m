classdef VolumeConstraintRhoD < Volume_constraint

    methods (Access = public)
        function computeFunctionAndGradient(obj)
            obj.computeFunctionAndGradient@Volume_constraint;
        end
    end
end