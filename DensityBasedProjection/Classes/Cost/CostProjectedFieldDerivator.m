classdef CostProjectedFieldDerivator < handle 
    %ProjectedCostDerivator
    properties (Access = public)
        derivedCost
    end
    properties (Access = private)
        structure
        mesh 
        displacement 
        projectedField
    end
    methods (Access = public)
        function obj = CostProjectedFieldDerivator(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.computeCost();
            obj.penalizeDerivatedCost();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.mesh.conectivityMatrixMat = cParams.conectivityMatrixMat;
            obj.mesh.elementNumberX = cParams.elementNumberX;
            obj.mesh.elementNumberY = cParams.elementNumberY;
            obj.structure = cParams.structure;
            obj.displacement = cParams.displacement; 
            obj.projectedField = cParams.projectedField;
            
        end 
        function computeCost(obj) 
           obj.derivedCost  = reshape(sum((obj.displacement(obj.mesh.conectivityMatrixMat)*obj.structure.elementalStiffnessMatrix).*obj.displacement(obj.mesh.conectivityMatrixMat),2),obj.mesh.elementNumberY,obj.mesh.elementNumberX);
        end 
        function penalizeDerivatedCost(obj)
            s.elasticModuleMinimun = obj.structure.elasticModuleMinimun;
            s.elasticModuleNeutral = obj.structure.elasticModuleNeutral;
            s.penalization = obj.structure.penalization;

            s.nonPenalizedVariable =  obj.derivedCost;
            s.projectedField = obj.projectedField ;
            B = DerivativePenalizer(s);
            B.penalize();
            obj.derivedCost = B.penalizedDerivative;
        end
         
    end
end