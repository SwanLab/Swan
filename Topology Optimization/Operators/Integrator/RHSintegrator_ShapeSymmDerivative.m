classdef RHSintegrator_ShapeSymmDerivative < RHSintegrator
    
    methods (Access = public)
        
        function obj = RHSintegrator_ShapeSymmDerivative(cParams)
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
            fG = fun.evaluate(); %% SIGMA
            dV = obj.mesh.computeDvolume(obj.quadrature);
            dNdx = test.computeCartesianDerivatives(obj.quadrature);
            nDim  = size(dNdx,1);
            nNode = size(dNdx,2);
            nElem = size(dNdx,3);
            nGaus = size(dNdx,4);

            BComp = obj.createBComputer(test,dNdx);
            rhsC = zeros(nNode*nDim,nElem);
            for igaus = 1:nGaus
                    fGI = squeezeParticular(fG(:,igaus,:),2);
                    fdv = fGI.*dV(igaus,:);
                    fdv = permute(fdv,[1 3 2]);      % (c x a x Z)
                    B = BComp.compute(igaus);        % (c x b x Z)

                    fp = permute(fdv,[1,2,4,3]);     % (c x a x 1 x Z)
                    Bp = permute(B,[1,4,2,3]);       % (c x 1 x b x Z)
                    intI = fp .* Bp;                 % (c x a x b x Z)
                    intI = sum(intI,1);              % (1 x a x b x Z)
                    intI = permute(intI,[2,3,4,1]);
                    rhsC = rhsC + squeezeParticular(intI,1);
            end
        end

        function f = assembleIntegrand(obj, rhsElem, test)
            integrand = rhsElem';
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