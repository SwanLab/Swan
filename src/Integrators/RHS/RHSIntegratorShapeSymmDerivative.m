classdef RHSIntegratorShapeSymmDerivative < RHSIntegrator
    
    methods (Access = public)
        
        function obj = RHSIntegratorShapeSymmDerivative(cParams)
            obj.init(cParams)
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
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end
        
        function rhsC = computeElementalRHS(obj, fun, test)
            xV    = obj.quadrature.posgp;
            fG    = fun.evaluate(obj.quadrature.posgp);
            symN  = ShapeDerSym(test).evaluate(xV);
            dV    = obj.mesh.computeDvolume(obj.quadrature);
            nDimf  = test.ndimf;
            nnodeE = obj.mesh.nnodeElem;
            nElem  = obj.mesh.nelem;
            ndofE  = nnodeE*nDimf;
            rhsC   = zeros(ndofE,nElem);
            for i=1:ndofE
                symTest = squeezeParticular(symN(:,:,i,:,:),3);
                fint    = pagetensorprod(fG,symTest,[1 2],[1 2],2,2);
                fint    = fint.*dV;
                rhsC(i,:) = rhsC(i,:) + sum(fint,1);
            end
        end

        function f = assembleIntegrand(obj, rhsElem, test)
            integrand = pagetranspose(rhsElem);
            connec = test.getDofConnec();
            nDofs = max(max(connec));
            nDofElem  = size(connec,2);
            f = zeros(nDofs,1);
            for idof = 1:nDofElem
                int = integrand(:,idof);
                con = connec(:,idof);
                f = f + accumarray(con,int,[nDofs,1],@sum,0);
            end
        end
        
    end
    
end