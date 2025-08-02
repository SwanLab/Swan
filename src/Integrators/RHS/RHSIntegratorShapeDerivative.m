classdef RHSIntegratorShapeDerivative < RHSIntegrator
    
    properties (Access = private)
        test
    end

    methods (Access = public)

        function obj = RHSIntegratorShapeDerivative(cParams)
            obj.init(cParams);
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
        end

        function rhs = compute(obj,fun)
            rhsElem = obj.computeElementalRHS(fun);
            rhs = obj.assembleIntegrand(rhsElem);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh            = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
            obj.test            = cParams.test;
        end
        
        function rhs = computeElementalRHS(obj,fun)
            xV     = obj.quadrature.posgp;
            dN     = ShapeDer(obj.test).evaluate(xV);
            fG     = fun.evaluate(xV);
            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.test.ndimf;
            ndofE  = nnodeE*ndim;
            nElem  = obj.mesh.nelem;
            rhs    = zeros(ndofE,nElem);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            for i = 1:ndofE
                dTest    = squeezeParticular(dN(:,:,:,:,i),1);
                intI     = pagetensorprod(fG,dTest,[1],[1],1,1);
                fI       = intI.*dV;
                rhs(i,:) = rhs(i,:) + sum(fI,1);
            end
            rhs = transpose(rhs);
        end

        function f = assembleIntegrand(obj,rhsElem)
            integrand = rhsElem;
            connec = obj.test.getDofConnec();
            nDofs = max(max(connec));
            nNode  = size(connec,2);
            f = zeros(nDofs,1);
            for inode = 1:nNode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[nDofs,1],@sum,0);
            end         
        end

    end

end