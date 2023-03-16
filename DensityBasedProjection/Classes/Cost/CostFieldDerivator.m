classdef CostFieldDerivator < handle
    properties (Access = public)
        derivedCost
    end
    properties (Access = private)
        mesh
        displacement
        structure
        projectedField
        projectorParameters
        derivedProjectedField
        filteredField
        filterParameters
        derivedFilteredField
        cost
    end
    methods (Access = public)
        function obj = CostFieldDerivator(cParams)
            obj.inputData(cParams);
        end

        function compute(obj)
            %Get the cost derivated respective the proyectedField
            obj.deriveCostRespectProjectedField();
            %Derivate the filtered field by the field
            obj.deriveFilteredFieldRespectedField();
            % Calculate the cost derivated by the field
            obj.derive();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.displacement =cParams.displacement;
            obj.structure =cParams.structure;
            obj.projectedField =cParams.projectedField;
            obj.projectorParameters =cParams.projectorParameters ;
            obj.filteredField =cParams.filteredField ;
            obj.filterParameters =cParams.filterParameters ;
            obj.cost =cParams.cost ;
            obj.derivedProjectedField = cParams.derivedProjectedField; 
            obj.derivedCost = zeros(obj.mesh.elementNumberX,obj.mesh.elementNumberY);

        end
        function derive(obj)
            obj.derivedCost(:) = obj.derivedFilteredField*(obj.cost.projectedDerived(:).*obj.derivedProjectedField(:));
        end
        function deriveCostRespectProjectedField(obj)
            s.elementNumberX = obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;
            s.structure = obj.structure;
            s.conectivityMatrixMat = obj.mesh.conectivityMatrixMat;
            s.projectedField = obj.projectedField;

            s.displacement =obj.displacement;
            B = CostProjectedFieldDerivator(s);
            B.compute();
            obj.cost.projectedDerived = B.derivedCost;
        end
        function deriveFilteredFieldRespectedField(obj)
            s.H = obj.filterParameters.H;
            s.Hs = obj.filterParameters.Hs;
            B = FilteredFieldFieldDerivator(s);
            B.compute();
            obj.derivedFilteredField = B.derivedFilteredField;
        end
    end
end