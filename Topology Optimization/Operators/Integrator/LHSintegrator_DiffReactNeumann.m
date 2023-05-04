classdef LHSintegrator_DiffReactNeumann < handle

    properties (GetAccess = public, SetAccess = private)
        M
        K
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = LHSintegrator_DiffReactNeumann(cParams)
            obj.mesh  = cParams.mesh;
            obj.computeStiffnessMatrix(cParams);
            obj.computeMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            LHS = epsilon^2*obj.K + obj.M;
        end

    end

    methods (Access = private)

        function computeStiffnessMatrix(obj,cParams)
            s      = cParams;  % For anisotropic stiffness
            s.type = cParams.stiffType;
            s.mesh = obj.mesh;
            s.fun  = P1Function.create(obj.mesh,1);
            LHS    = LHSintegrator.create(s);
            obj.K  = LHS.compute();
        end

        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATICMASS';
            LHS     = LHSintegrator.create(s);
            obj.M   = LHS.compute();
        end

    end

end