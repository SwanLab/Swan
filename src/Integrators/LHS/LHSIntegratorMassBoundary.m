classdef LHSIntegratorMassBoundary < handle

    properties (Access = private)
        mesh
        trial
        test
    end

    methods (Access = public)
        function obj = LHSIntegratorMassBoundary(cParams)
            obj.mesh  = cParams.mesh;
            obj.trial = cParams.trial;
            obj.test  = cParams.test;
        end

        function LHS = compute(obj)
            [bMesh, l2g] = obj.mesh.createSingleBoundaryMesh();
            ndof         = obj.trial.nDofs;
            LHS          = sparse(ndof,ndof);
            Mloc         = obj.computeLocalMassMatrix(bMesh);
            LHS(l2g,l2g) = Mloc;
        end
    end

    methods (Access = private)
        function M = computeLocalMassMatrix(obj,m)
            s.type  = 'MassMatrix';
            s.mesh  = m;
            s.test  = LagrangianFunction.create(m,obj.test.ndimf,obj.test.order);
            s.trial = LagrangianFunction.create(m,obj.trial.ndimf,obj.trial.order);
            lhs     = LHSIntegrator.create(s);
            M       = lhs.compute();
        end
    end
end