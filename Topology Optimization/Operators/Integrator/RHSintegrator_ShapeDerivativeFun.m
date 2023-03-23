classdef RHSintegrator_ShapeDerivativeFun < handle

    properties (Access = private)
        npnod
        mesh
        fNodal
        quadratureOrder
        quadrature

        fun
    end

    methods (Access = public)

        function obj = RHSintegrator_ShapeDerivativeFun(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function rhs = compute(obj)
            rhsElem = obj.computeElementalRHS();
            rhs = obj.assembleIntegrand(rhsElem);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.npnod        = cParams.npnod;
            obj.quadratureOrder = cParams.quadratureOrder;
            obj.fun      = cParams.fun;
        end

        function createQuadrature(obj)
            q = Quadrature.create(obj.mesh, obj.quadratureOrder);
            obj.quadrature = q;
        end
        
        function rhsC = computeElementalRHS(obj)
            fG    = obj.fun.evaluate(obj.quadrature.posgp);
            dV    = obj.mesh.computeDvolume(obj.quadrature);
            dNdx  = obj.fun.computeCartesianDerivatives();
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
            ndofs = obj.npnod;
            connec = obj.mesh.connec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end

    end

end