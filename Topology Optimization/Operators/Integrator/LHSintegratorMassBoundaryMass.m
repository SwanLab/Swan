classdef LHSintegratorMassBoundaryMass < handle

    properties (GetAccess = public, SetAccess = private)
        test, trial
        M
        Mr
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        function obj = LHSintegratorMassBoundaryMass(cParams)
            obj.init(cParams);
            obj.computeMassMatrix();
            obj.computeBoundaryMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            LHS = obj.M + epsilon*obj.Mr;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.test  = cParams.trial;
            obj.trial = cParams.trial;
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