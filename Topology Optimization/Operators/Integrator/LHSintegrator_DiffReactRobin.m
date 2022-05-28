classdef LHSintegrator_DiffReactRobin < LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        M
        K
        Mr
    end

    properties (Access = private)
        field
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
            g.mesh               = obj.mesh;
            g.ndimf              = 1;
            g.interpolationOrder = 'LINEAR';
            f = Field(g);
            s.type = 'StiffnessMatrix';
            s.mesh  = obj.mesh;
            s.field = f;
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();
        end

        function computeMassMatrix(obj)
            g.mesh               = obj.mesh;
            g.ndimf              = 1;
            g.interpolationOrder = 'LINEAR';
            g.quadratureOrder    = 'QUADRATICMASS';
            f = Field(g);
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = f;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end

        function computeBoundaryMassMatrix(obj)
            g.mesh               = obj.mesh;
            g.ndimf              = 1;
            g.interpolationOrder = 'LINEAR';
            g.quadratureOrder    = 'QUADRATICMASS';
            f = Field(g);
            s.type         = 'BoundaryMassMatrix';
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.globalConnec = [];
            s.field = f;
            LHS = LHSintegrator.create(s);
            obj.Mr = LHS.compute();
        end

    end

end