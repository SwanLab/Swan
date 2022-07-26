classdef LHSintegrator_DiffReactNeumann < LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        M
        K
    end

    properties (Access = private)
        field
    end

    methods (Access = public)

        function obj = LHSintegrator_DiffReactNeumann(cParams)
            obj.mesh  = cParams.mesh;
            obj.field = cParams.field;
            obj.computeStiffnessMatrix(cParams);
            obj.computeMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            LHS = epsilon^2*obj.K + obj.M;
        end

    end

    methods (Access = private)

        function computeStiffnessMatrix(obj,cParams)
            s              = cParams;
            s.globalConnec = obj.mesh.connec;
            s.type         = cParams.stiffType;
            s.mesh         = obj.mesh;
            s.field        = obj.field;
            LHS            = LHSintegrator.create(s);
            obj.K          = LHS.compute();
        end

        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            LHS     = LHSintegrator.create(s);
            obj.M   = LHS.compute();
        end

    end

end