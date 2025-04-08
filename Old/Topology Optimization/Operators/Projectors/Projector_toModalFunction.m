classdef Projector_toModalFunction < Projector

    properties (Access = private)
        basis
        functionType
    end

    methods (Access = public)

        function obj = Projector_toModalFunction(cParams)
            obj.init(cParams);
            obj.basis = cParams.basis;
            obj.functionType = cParams.functionType;
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            s.basis = obj.basis;
            s.functionType = obj.functionType;
            xFun = ModalFunction(s);
        end

    end

    methods (Access = private)

        function LHS = computeLHS(obj)
            mesh         = obj.mesh;
            basis        = obj.basis;
            functionType = obj.functionType;

            test  = ModalFunction.create(mesh,basis,functionType);
            quad  = obj.createRHSQuadrature(test);
            xV    = quad.posgp;
            dV    = obj.mesh.computeDvolume(quad);
            ngaus = quad.ngaus;

            basisTest  = test.evaluateBasisFunctions(xV);
            basisTrial = basisTest;
            nbasis     = test.nbasis;

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
            basis    = obj.basis;
            functionType = obj.functionType;

            test  = ModalFunction.create(mesh,basis,functionType);

            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);

            basisTest  = test.evaluateBasisFunctions(xV);
            nbasis     = test.nbasis;

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