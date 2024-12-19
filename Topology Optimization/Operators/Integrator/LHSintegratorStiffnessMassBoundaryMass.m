classdef LHSintegratorStiffnessMassBoundaryMass < handle

    properties (GetAccess = public, SetAccess = private)
        test, trial
        M
        K
        Mr
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = LHSintegratorStiffnessMassBoundaryMass(cParams)
            obj.init(cParams);
            obj.computeStiffnessMatrix(cParams);
            obj.computeMassMatrix();
            obj.computeBoundaryMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            lump = sum(obj.M + epsilon*obj.Mr,2);
            LHS  = epsilon^2*obj.K + diag(lump);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.test  = cParams.trial;
            obj.trial = cParams.trial;
        end

        function computeStiffnessMatrix(obj,cParams)
            s       = cParams;
            s.test  = obj.test;
            s.trial = obj.trial;
            s.type  = cParams.stiffType;
            s.mesh  = obj.mesh;
            LHS     = LHSintegrator.create(s);
            obj.K   = LHS.compute();
        end

        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.test;
            s.trial = obj.trial;
            s.quadratureOrder = 2;
            LHS     = LHSintegrator.create(s);
            obj.M   = LHS.compute();
        end

        function computeBoundaryMassMatrix(obj)
            s.type  = 'BoundaryMassMatrix';
            s.mesh  = obj.mesh;
            s.quadratureOrder = 2;
            LHS     = LHSintegrator.create(s);
            obj.Mr  = LHS.compute();
        end

    end

end