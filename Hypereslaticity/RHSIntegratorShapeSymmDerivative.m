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
            xV     = obj.quadrature.posgp;
            dSymN  = ShapeDerSym(test);
            symN   = dSymN.evaluate(xV);
            fG     = fun.evaluate(xV);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            nnodeE = obj.mesh.nnodeElem;
            ndim   = test.ndimf;
            ndofE  = nnodeE*ndim;
            nElem  = obj.mesh.nelem;
            rhsC   = zeros(ndofE,nElem);
            for i = 1:ndofE
                    symTest   = squeezeParticular(symN(:,:,i,:,:),3);
                    dfint     = pagetensorprod(fG,symTest,[1 2],[1 2],2,2);
                    dfint     = dfint.*dV;
                    rhsC(i,:) = rhsC(i,:) + sum(dfint,1);
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

        function BComp = createBComputer(obj, fun, dNdx)
            s.fun = fun;
            s.dNdx = dNdx;
            BComp = BMatrixComputer(s);
        end
    end
    
end