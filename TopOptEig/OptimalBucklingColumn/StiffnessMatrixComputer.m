classdef StiffnessMatrixComputer < handle

    properties (Access = private)
        LHS
        stiffnessMatrix
    end

    properties (Access = private)
        dim
        mesh
        youngModulus
        inertiaMoment
        freeNodes
    end
    
    methods (Access = public)

        function obj = StiffnessMatrixComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.createStiffnessMatrix();
            obj.computeStiffnessMatrix()
        end
       
        function Kfree = provideFreeStiffnessMatrix(obj)
            free = obj.freeNodes;
            K = obj.stiffnessMatrix;
            Kfree  = K(free,free);
        end

    end
    
    methods (Access = private)
        
        function obj = init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.dim           = cParams.dim;
            obj.freeNodes     = cParams.freeNodes;
        end

        function createStiffnessMatrix(obj)
            s.type = 'StiffnessMatrixColumn';
            s.dim = obj.dim;
            s.mesh = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            obj.LHS = LHSintegrator.create(s);
        end

        function computeStiffnessMatrix(obj)
            obj.stiffnessMatrix = obj.LHS.compute();
        end
 
    end
    
end