classdef Projector_toRigidBody < Projector

    properties (Access = private)
        refPoint
    end

    methods (Access = public)

        function obj = Projector_toRigidBody(cParams)
            obj.init(cParams);
            obj.refPoint = cParams.refPoint;
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fvalues = xProj;
            s.refPoint = obj.refPoint;
            xFun = RigidBodyFunction(s);
        end

    end

    methods (Access = private)

        function LHS = computeLHS(obj)
            mesh     = obj.mesh;
            refPoint = obj.refPoint;
            test  = RigidBodyFunction.create(mesh,refPoint);
            quad  = obj.createRHSQuadrature(test);
            xV    = quad.posgp;
            dV    = obj.mesh.computeDvolume(quad);
            ngaus = quad.ngaus;
            basisTest  = test.computeBasisFunction(xV);
            basisTrial = basisTest;
            nbasis     = size(basisTest,2);
            LHS   = zeros(nbasis);
            nFlds = test.ndimf;
            for i = 1:nbasis
                for j = 1:nbasis
                    for iField = 1:nFlds
                        basisProd = squeeze(basisTest{i}(iField,:,:).*basisTrial{j}(iField,:,:));
                        LHS(i,j)  = LHS(i,j) + sum(basisProd.*dV,'all');
                    end
                end
            end
        end

        function RHS = computeRHS(obj,fun)
            mesh     = obj.mesh;
            refPoint = obj.refPoint;

            test  = RigidBodyFunction.create(mesh,refPoint);

            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);

            basisTest  = test.computeBasisFunction(xV);
            nbasis     = size(basisTest,2);

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;

            fGaus = fun.evaluate(xV);
            f     = zeros(nbasis,1);

            for i = 1:nbasis
                for iField = 1:nFlds
                    int = squeeze(basisTest{i}(iField,:,:).*fGaus(iField,:,:));
                    f(i) = f(i) + sum(int.*dV,'all');
                end
            end

            RHS = f;
        end

        function q = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end

    end

end