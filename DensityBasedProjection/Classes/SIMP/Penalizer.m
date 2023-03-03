classdef Penalizer < SIMPcomputer
    properties (Access = public)
        penalizedVariable
    end
    methods (Access = public)
        function penalize(obj)
            obj.penalizedVariable = obj.nonPenalizedVariable(:)*(obj.elasticModuleMinimun + obj.projectedField(:)'.^obj.penalization*(obj.elasticModuleNeutral-obj.elasticModuleMinimun));
        end
    end
end