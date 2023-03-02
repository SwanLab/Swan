classdef LHSintegrator_DiffReactRobin < LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        M
        K
        Mr
    end

    methods (Access = public)

        function obj = LHSintegrator_DiffReactRobin(cParams)
            obj.mesh  = cParams.mesh;
            obj.computeStiffnessMatrix(cParams);
            obj.computeMassMatrix();
            obj.computeBoundaryMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            LHS = epsilon^2*obj.K + obj.M + epsilon*obj.Mr;
        end

    end

    methods (Access = private)

        function computeStiffnessMatrix(obj,cParams)
            g.mesh               = obj.mesh;
            g.ndimf              = 1;
            g.interpolationOrder = 'LINEAR';
            g.quadratureOrder    = 'LINEAR';
            f = Field(g);
            s              = cParams;
            s.globalConnec = obj.mesh.connec;
            s.type         = cParams.stiffType;
            s.mesh         = obj.mesh;
            s.field        = f;
            LHS            = LHSintegrator.create(s);
            obj.K          = LHS.compute();
        end

        function computeMassMatrix(obj)
            g.mesh    = obj.mesh;
            g.fValues = zeros(size(obj.mesh.coord,1),1);
            f = P1Function(g);
            s.type = 'MassMatrixFun';
            s.mesh = obj.mesh;
            s.fun  = f;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS     = LHSintegrator.create(s);
            obj.M   = LHS.compute();
        end

        function computeBoundaryMassMatrix(obj)
            g.mesh    = obj.mesh;
            g.fValues = zeros(size(obj.mesh.coord,1),1);
            f = P1Function(g);
            s.type  = 'BoundaryMassMatrix';
            s.mesh  = obj.mesh;
            s.fun   = f;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS     = LHSintegrator.create(s);
            obj.Mr  = LHS.compute();
        end

    end

end