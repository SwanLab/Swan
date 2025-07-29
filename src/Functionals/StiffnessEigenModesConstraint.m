classdef StiffnessEigenModesConstraint < handle

    properties (Access = private)
        mesh
        eigenModesFunctional
        designVariable
        targetEigenValue
        filteredDesignVariable
        iter
        filter
        value0
        material
    end
    
    methods (Access = public)
        function obj = StiffnessEigenModesConstraint(cParams)
            obj.init(cParams);
            eigen = StiffnessEigenModesDisplacementComputer(cParams);
            s = cParams;
            s.eigenModes = eigen;
            obj.eigenModesFunctional = MinimumEigenValueFunctional(s);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.targetEigenValue = cParams.targetEigenValue;
            obj.designVariable    = cParams.designVariable;
            obj.iter              = 0.0;
            obj.material          = cParams.material;
        end
    end
end

