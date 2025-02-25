classdef BendingMatrixComputer < handle
    
    properties (Access = public)
        elementalBendingMatrix
    end

    properties (Access = private)
        bendingMatrix
        LHS
    end

    properties (Access = private)
        dim
        mesh
        youngModulus
        inertiaMoment
        designVariable
        freeNodes
    end
    
    methods (Access = public)
        
        function obj = BendingMatrixComputer(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.createBendingMatrix();
            obj.computeBendingMatrix();
        end

        function Bfree = provideFreeBendingMatrix(obj)
            free = obj.freeNodes;
            B = obj.bendingMatrix;
            Bfree  = B(free,free);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.dim            = cParams.dim;
            obj.mesh           = cParams.mesh;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.designVariable = cParams.designVariable;
            obj.freeNodes      = cParams.freeNodes;
        end

        function createBendingMatrix(obj)
            s.type         = 'BendingMatrix';
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.inertiaMoment  = obj.inertiaMoment;
            s.youngModulus   = obj.youngModulus;
            s.designVariable = obj.designVariable;
            obj.LHS = LHSintegrator.create(s);
        end

        function computeBendingMatrix(obj)
            obj.bendingMatrix = obj.LHS.compute();
            obj.elementalBendingMatrix = obj.LHS.elementalBendingMatrix;
        end

    end
end