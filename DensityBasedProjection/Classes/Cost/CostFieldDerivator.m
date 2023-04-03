classdef CostFieldDerivator < handle
    properties (Access = public)
        derivedCost
    end
    properties (Access = private)
        mesh
        displacement
        structure
        projectorParameters
        designFields
        cost
        derivedCostProjected
    end
    methods (Access = public)
        function obj = CostFieldDerivator(cParams)
            obj.inputData(cParams);
        end

        function compute(obj)
            %Get the cost derivated respective the proyectedField
            obj.deriveCostRespectProjectedField();
            %Penalize the the cost derivated respective the proyectedField
            obj.penalizeDerivatedCost();
            % Calculate the cost derivated by the field
            obj.derive();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.displacement =cParams.displacement;
            obj.structure =cParams.structure;
            obj.designFields = cParams.designFields;
            obj.cost =cParams.cost ;
            obj.derivedCost = zeros(obj.mesh.elementNumberX,obj.mesh.elementNumberY);

        end
        function derive(obj)
            obj.derivedCost(:) = obj.designFields.derivedFilteredField*(obj.derivedCostProjected(:).*obj.designFields.derivedProjectedField(:));
        end
        function deriveCostRespectProjectedField(obj)
            s.elementNumberX = obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;
            s.structure = obj.structure;
            s.conectivityMatrixMat = obj.mesh.conectivityMatrixMat;
            s.projectedField = obj.designFields.projectedField;

            s.displacement =obj.displacement;
            B = CostProjectedFieldDerivator(s);
            B.compute();
            obj.derivedCostProjected = B.derivedCost;
        end
        function penalizeDerivatedCost(obj)
            s.elasticModuleMinimun = obj.structure.elasticModuleMinimun;
            s.elasticModuleNeutral = obj.structure.elasticModuleNeutral;
            s.penalization = obj.structure.penalization;
            s.nonPenalizedVariable =  obj.derivedCostProjected;
            s.projectedField = obj.designFields.projectedField ;
            B = DerivativePenalizer(s);
            B.penalize();
            obj.derivedCostProjected = B.penalizedDerivative;
        end
    end
end