classdef LHSintegrator_DiffReactNeumann < LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        M
        K
    end

    methods (Access = public)

        function obj = LHSintegrator_DiffReactNeumann(cParams)
            obj.init(cParams);
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            LHS = epsilon^2*obj.K + obj.M;
        end

    end

    methods (Access = private)
    
        function computeStiffnessMatrix(obj)
            s.type = 'StiffnessMatrix';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();
        end
        
        function computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end


    end

end