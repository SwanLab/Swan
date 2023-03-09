classdef SIMPcomputer < handle

    properties (Access = public)
        elasticModuleMinimun
        elasticModuleNeutral
        projectedField
        penalization
        nonPenalizedVariable
    end

    methods (Access = public)
        function obj = SIMPcomputer(cParams)
            obj.inputData(cParams); 
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.elasticModuleMinimun = cParams.elasticModuleMinimun;
            obj.elasticModuleNeutral = cParams.elasticModuleNeutral;
            obj.nonPenalizedVariable = cParams.nonPenalizedVariable;
            obj.projectedField = cParams.projectedField;
            obj.penalization = cParams.penalization;
        end
    end
end