classdef ComplianceConstraintThreeField < ShFunc_Compliance

    methods (Access = protected)
        function normalizeFunction(obj)
            obj.normalizeFunction@ShFunc_Compliance;
            obj.value = obj.value/2 - 1;
        end
    end
end