classdef RHSintegrator_ShapeDivergence < handle

    properties (Access = private)
        mesh
        quadratureOrder
        quadrature
        testF
    end

    methods (Access = public)

        function obj = RHSintegrator_ShapeDivergence(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function rhsFun = compute(obj, fun)
            rhsElem = obj.computeElementalRHS(fun);
            for iDim = 1:obj.mesh.ndim
            rhs(:,iDim) = obj.assembleIntegrand(rhsElem(:,:,iDim));
            end
            s.fValues = rhs;
            s.mesh    = obj.mesh;
            rhsFun = P1Function(s);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh            = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
            s.mesh              = obj.mesh;
            s.fValues           = [];
            obj.testF           = P1Function(s);
        end

        function createQuadrature(obj)
            q = Quadrature.create(obj.mesh, obj.quadratureOrder);
            obj.quadrature = q;
        end
        
     

       function rhsC = computeElementalRHS(obj, fun)
            fG    = fun.evaluate(obj.quadrature.posgp);
            dV    = obj.mesh.computeDvolume(obj.quadrature);
            dNdx = obj.testF.computeCartesianDerivatives(obj.quadrature);
            nDim  = size(dNdx,1);
            nNode = size(dNdx,2);
            nElem = size(dNdx,3);
            nGaus = size(dNdx,4);
            int = zeros(nNode,nElem,nDim);
            for idime = 1:nDim            
               for igaus = 1:nGaus
                    for inode = 1:nNode
                        fI     = squeezeParticular(fG(1,igaus,:),1);
                        fdV    = fI.*dV(igaus,:);
                        dShape = squeeze(dNdx(idime,inode,:,igaus))';
                        intI = dShape.*fdV;
                        int(inode,:,idime) = int(inode,:,idime) + intI;
                    end
                end
            end
            rhsC = permute(int,[2 1 3]);
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