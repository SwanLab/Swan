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
%         function computeSIMPelasticModule(obj)
%             obj.SIMPelasticModule = obj.elasticModuleMinimun + obj.xPhys(:)'.^obj.penalization*(obj.elasticModuleNeutral-obj.elasticModuleMinimun);
%         end
%         function computeDerivativeSIMPelasticModule(obj)
%             obj.derivativeSIMPelasticModule = -obj.penalization*(obj.elasticModuleNeutral-obj.elasticModuleMinimun)*xPhysE.^(obj.penalization-1).*cE1a;
%         end
%         function compute(obj)
%             obj.penalize();
%         end 
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