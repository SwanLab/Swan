classdef LHSIntegratorStiffnessMass < handle

    properties (GetAccess = public, SetAccess = private)
        test, trial
        M
        K
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = LHSIntegratorStiffnessMass(cParams)
            obj.init(cParams);
            obj.computeStiffnessMatrix(cParams);
            obj.computeMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            LHS = epsilon^2*obj.K + obj.M;
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
            s.type  = cParams.stiffType;
            s.mesh  = obj.mesh;
            s.test  = obj.test;
            s.trial = obj.trial;
            s.quadratureOrder = 2;
            LHS     = LHSIntegrator.create(s);
            obj.K   = LHS.compute();
        end

        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.test;
            s.trial = obj.trial;
            s.quadratureOrder = 2;
            LHS     = LHSIntegrator.create(s);
            obj.M   = LHS.compute();
        end

    end

end