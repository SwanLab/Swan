classdef FEMcomputer < handle
    properties (Access = public)
        displacement
        cost
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
            obj.computeStifnessMatrix
            obj.computeForces
            obj.computeDisplacement;
            obj.computeCost;

        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.structure =cParams.structure;
            obj.projectedField = cParams.projectedField;
            obj.displacement = zeros(2*(obj.mesh.elementNumberY+1)*(obj.mesh.elementNumberX+1),1);
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
            s.neumanCondition = obj.mesh.neumanCondition;
            s.output = obj.mesh.output;
            s.elementNumberX = obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;
            B = ForceComputer(s);
            B.compute();
            obj.force = B.force;
        end
        function computeDisplacement(obj)
            s.force = obj.force;
            s.globalStifnessMatrix =  obj.globalStifnessMatrix;
            s.freeDegress = obj.mesh.freeDegress;
            s.elementNumberX = obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;
            B = DisplacementComputer(s);
            B.compute();
            obj.displacement = B.displacement;
        end
        function computeCost(obj)
            s.force = obj.force;
            s.displacement =obj.displacement;
            B = CostComputer(s);
            B.compute;
            obj.cost = B.cost;
        end
    end
end