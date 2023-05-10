classdef ComplianceConstraintThreeField < ShFunc_Compliance

    methods (Access = protected)
        function computeFunctionValue(obj)
            obj.computeFunctionValue@ShFunc_Compliance;
            % ... obj.value = obj.value/target - 1
        end
    end
end