classdef FEMcomputer < handle
    properties (Access = public)
        displacement
        force
    end
    properties (Access = private)
        mesh
        structure
        projectedField
        globalStifnessMatrix
    end

    methods (Access = public)
        function obj = FEMcomputer(cParams)
            obj.inputData(cParams)
        end
        function compute(obj)
            obj.penalizeElasticModule
            obj.computeStifnessMatrix
            obj.computeForces
            obj.computeDisplacement;
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.structure =cParams.structure;
            obj.projectedField = cParams.projectedField;
            obj.displacement = zeros(2*(obj.mesh.elementNumberY+1)*(obj.mesh.elementNumberX+1),1);
        end
        function penalizeElasticModule(obj)
            s.elasticModuleMinimun = obj.structure.elasticModuleMinimun;
            s.elasticModuleNeutral = obj.structure.elasticModuleNeutral;
            s.projectedField = obj.projectedField;
            s.penalization = obj.structure.penalization;
            s.nonPenalizedVariable = obj.structure.elementalStiffnessMatrix;
            B = Penalizer(s);
            B.penalize();
            obj.structure.elementalStiffnessMatrix = B.penalizedVariable;
        end
        function computeStifnessMatrix(obj)
            s.structure = obj.structure;
            s.mesh = obj.mesh;
            s.projectedField            =  obj.projectedField;
            B = StifnessMatrixComputer(s);
            B.compute();
            obj.globalStifnessMatrix = B.globalStifnessMatrix;
        end

        function computeForces(obj)
            s.allDegrees = obj.mesh.degress.all;
            s.neumanCondition = obj.mesh.degress.forceDOFs(:,2); 
            s.output =  obj.mesh.degress.forceDOFs(:,1);
            s.elementNumberX = obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;
            B = ForceComputer(s);
            B.compute();
            obj.force = B.force;
        end
        function computeDisplacement(obj)
            s.force = obj.force;
            s.globalStifnessMatrix =  obj.globalStifnessMatrix;
            s.freeDegress = obj.mesh.degress.free;
            s.elementNumberX = obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;
            B = DisplacementComputer(s);
            B.compute();
            obj.displacement = B.displacement;
        end
    end
end