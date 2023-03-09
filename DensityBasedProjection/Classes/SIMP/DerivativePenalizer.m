classdef DerivativePenalizer < SIMPcomputer
    properties (Access = public)
        penalizedDerivative
    end
    properties (Access = private)
        derivativeCost
    end
    methods (Access = public)
        function penalize(obj)
            obj.penalizedDerivative = -obj.penalization*(obj.elasticModuleNeutral-obj.elasticModuleMinimun)*obj.projectedField.^(obj.penalization-1).*obj.nonPenalizedVariable;
        end
    end
end
