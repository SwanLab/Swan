classdef RHSintegrator_ShapeDerivative < handle

    properties (Access = private)
        mesh
        quadratureOrder
        quadrature
    end

    methods (Access = public)

        function obj = RHSintegrator_ShapeDerivative(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function rhs = compute(obj, fun)
            rhsElem = obj.computeElementalRHS(fun);
            rhs = obj.assembleIntegrand(rhsElem);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function createQuadrature(obj)
            q = Quadrature.create(obj.mesh, obj.quadratureOrder);
            obj.quadrature = q;
        end
        
        function rhsC = computeElementalRHS(obj, fun)
            fG    = fun.evaluate(obj.quadrature.posgp);
            dV    = obj.mesh.computeDvolume(obj.quadrature);
            dNdx  = fun.computeCartesianDerivatives(obj.quadrature);
            nDim  = size(dNdx,1);
            nNode = size(dNdx,2);
            nElem = size(dNdx,3);
            nGaus = size(dNdx,4);
            int = zeros(nNode,nElem);
            for igaus = 1:nGaus
                for idime = 1:nDim
                    for inode = 1:nNode
                        fI     = squeezeParticular(fG(idime,igaus,:),1);
                        fdV    = fI.*dV(igaus,:);
                        dShape = squeeze(dNdx(idime,inode,:,igaus))';
                        intI = dShape.*fdV;
                        int(inode,:) = int(inode,:) + intI;
                    end
                end
            end
            rhsC = transpose(int);
        end

        function f = assembleIntegrand(obj,rhsElem)
            integrand = rhsElem;
            connec = obj.mesh.connec;
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