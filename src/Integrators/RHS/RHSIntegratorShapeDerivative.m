classdef RHSIntegratorShapeDerivative < RHSIntegrator

    properties (Access = private)
        test
        dNdx
    end

    methods (Access = public)

        function obj = RHSIntegratorShapeDerivative(cParams)
            obj.init(cParams);
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
            obj.computeShapeDerivatives();
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

        function computeShapeDerivatives(obj)
            xV       = obj.quadrature.posgp;
            dN       = obj.test.evaluateCartesianDerivatives(xV);
            obj.dNdx = permute(dN,[5 4 1 2 3]);
        end

        function rhsC = computeElementalRHS(obj,fun)
            xV = obj.quadrature.posgp;
            fG    = fun.evaluate(xV);
            dV    = obj.mesh.computeDvolume(obj.quadrature);
            nElem = size(obj.dNdx,2);
            nDim  = size(obj.dNdx,3);
            nNode = size(obj.dNdx,4);
            nGaus = size(obj.dNdx,5);
            int = zeros(nNode,nElem);
            for igaus = 1:nGaus
                for idime = 1:nDim
                    for inode = 1:nNode
                        fI     = squeezeParticular(fG(idime,igaus,:),[1 2]);
                        fdV    = fI'.*dV(igaus,:);
                        dShape = obj.dNdx(1,:,idime,inode,igaus);
                        intI = dShape.*fdV;
                        int(inode,:) = int(inode,:) + intI;
                    end
                end
            end
            rhsC = transpose(int);
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