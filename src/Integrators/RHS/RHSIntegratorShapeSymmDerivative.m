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
            fG = fun.evaluate(obj.quadrature.posgp);
            dV = obj.mesh.computeDvolume(obj.quadrature);
            xV     = obj.quadrature.posgp;
            % dNdx = test.evaluateCartesianDerivatives(obj.quadrature.posgp);
            dSymN  = ShapeDerSym(test);
            symN   = dSymN.evaluate(xV);
            nDimf  = test.ndimf;
            nnodeE = obj.mesh.nnodeElem;
            % nGaus = obj.
            nElem = obj.mesh.nelem;
            ndofE=nnodeE*nDimf;

            % BComp = obj.createBComputer(test,dNdx);
            rhsC = zeros(ndofE,nElem);
            % for igaus = 1:nGaus
            %         fGI = squeezeParticular(fG(:,igaus,:),2);
            %         fdv = fGI'.*dV(igaus,:);
            %         fdv = reshape(fdv,[1 size(fdv,1) nElem]);
            %         B = BComp.compute(igaus);
            %         intI = pagemtimes(fdv,B);
            %         rhsC = rhsC + squeezeParticular(intI,1);
            % end
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

        function BComp = createBComputer(obj, fun, dNdx)
            s.fun = fun;
            s.dNdx = dNdx;
            BComp = BMatrixComputer(s);
        end
    end
    
end