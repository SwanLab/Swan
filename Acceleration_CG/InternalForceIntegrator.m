classdef InternalForceIntegrator < RHSintegrator
    
    methods (Access = public)
        
        function obj = InternalForceIntegrator(cParams)
            obj.init(cParams)
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
        end
        
        function rhs = compute(obj, fun, test)
            rhsElem = obj.computeElementalRHS(fun,test);
            rhs     = obj.assembleIntegrand(rhsElem,test);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end
        
        function rhsC = computeElementalRHS(obj, fun, test)
            disp('-------')
            tic
            fG = fun.evaluate(obj.quadrature.posgp);
            toc
            tic
            dV = obj.mesh.computeDvolume(obj.quadrature);
            toc
            tic
            dNdx = test.evaluateCartesianDerivatives(obj.quadrature.posgp);
            toc
            tic
            nDim  = size(dNdx,1);
            nNode = size(dNdx,2);
            nGaus = size(dNdx,3);
            nElem = size(dNdx,4);

            BComp = obj.createBComputer(test,dNdx);
            rhsC = zeros(nNode*nDim,nElem);
            for igaus = 1:nGaus
                    fGI = squeezeParticular(fG(:,igaus,:),2);
                    fdv = fGI.*dV(igaus,:);
                    fdv = reshape(fdv,[1 size(fdv,1) nElem]);
                    % tic
                    B = BComp.compute(igaus);
                    % disp(strcat('B computation: ',string(toc),' s'))
                    intI = pagemtimes(fdv,B);
                    rhsC = rhsC + squeezeParticular(intI,1);
            end
            toc
            % disp(strcat('Time RHS assembly: ',string(toc),' s'))
        end

        function f = assembleIntegrand(obj, rhsElem, test)
            integrand = pagetranspose(rhsElem);
            connec = test.getConnec();
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