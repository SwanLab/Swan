classdef RHSintegrator_ShapeDerivative < RHSintegrator

    methods (Access = public)

        function obj = RHSintegrator_ShapeDerivative(cParams)
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
        
        function rhsC = computeElementalRHS(obj, fun, test)
            xV = obj.quadrature.posgp;
            fG    = fun.evaluate(xV);
            dNdx  = test.evaluateCartesianDerivatives(xV);
            dV    = obj.mesh.computeDvolume(obj.quadrature);
            nDim  = size(dNdx,1);
            nNode = size(dNdx,2);
            nGaus = size(dNdx,3);
            nElem = size(dNdx,4);
            int = zeros(nNode,nElem);
            for igaus = 1:nGaus
                for idime = 1:nDim
                    for inode = 1:nNode
                        fI     = squeezeParticular(fG(idime,igaus,:),1);
                        fdV    = fI.*dV(igaus,:);
                        dShape = squeeze(dNdx(idime,inode,igaus,:))';
                        intI = dShape.*fdV;
                        int(inode,:) = int(inode,:) + intI;
                    end
                end
            end
            rhsC = transpose(int);
        end

        function f = assembleIntegrand(obj,rhsElem,test)
            integrand = rhsElem;
            connec = test.getDofConnec();
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