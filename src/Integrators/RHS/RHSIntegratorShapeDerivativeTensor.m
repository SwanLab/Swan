classdef RHSIntegratorShapeDerivativeTensor < RHSIntegrator

    methods (Access = public)

        function obj = RHSIntegratorShapeDerivativeTensor(cParams)
            obj.init(cParams);
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
        end

        function rhs = compute(obj, fun, test)
            rhsElem = obj.computeElementalRHS(fun,test);
            rhs = obj.assembleIntegrand(rhsElem,test);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function rhs = computeElementalRHS(obj, fun, test)
            xV     = obj.quadrature.posgp;
            dN     = ShapeDer(test).evaluate(xV);
            fG     = fun.evaluate(xV);
            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.mesh.ndim;
            ndofE  = nnodeE*ndim;
            nElem  = obj.mesh.nelem;
            rhs    = zeros(ndofE,nElem);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            for i = 1:ndofE
                dTest = dN(:,:,:,:,i);
                intI  = pagetensorprod(fG,dTest,[1 2],[1 2],2,2);
                fI    = intI.*dV;
                rhs(i,:) = rhs(i,:) + sum(fI,1);
            end
            rhs = rhs';
        end

        function f = assembleIntegrand(obj,rhsElem,test)
            integrand = rhsElem;
            connec = test.getDofConnec();
            nDofs = max(max(connec));
            nDofE = size(connec,2);
            f = zeros(nDofs,1);
            for idof = 1:nDofE
                int = integrand(:,idof);
                con = connec(:,idof);
                f = f + accumarray(con,int,[nDofs,1],@sum,0);
            end         
        end

    end

end