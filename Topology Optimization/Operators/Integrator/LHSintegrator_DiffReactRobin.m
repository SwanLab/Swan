classdef LHSintegrator_DiffReactRobin < LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        M
        K
        Mr
    end

    methods (Access = public)

        function obj = LHSintegrator_DiffReactRobin(cParams)
            obj.init(cParams);
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix();
            obj.computeBoundaryMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            LHS = epsilon^2*obj.K + obj.M + epsilon*obj.Mr;
        end

    end

    methods (Access = private)
    
        function computeStiffnessMatrix(obj)
            s.type = 'StiffnessMatrix';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();
        end
        
        function computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end

        function computeBoundaryMassMatrix(obj)
            s.type         = 'BoundaryMassMatrix';
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.globalConnec = [];
            LHS = LHSintegrator.create(s);
            obj.Mr = LHS.compute();
        end

    end

end