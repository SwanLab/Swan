classdef RHSintegrator_ShapeSymmetricGradient < handle

    properties (Access = private)
        mesh
        quadratureOrder
        quadrature
        testF
    end

    methods (Access = public)

        function obj = RHSintegrator_ShapeSymmetricGradient(cParams)
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
     
       function rhsS = computeElementalRHS(obj, fun)
            fG   = fun.evaluate(obj.quadrature.posgp);
            rhsC = obj.computeElementalGradient(fG);
            fGT  = permute(fG,[2 1 3 4]);
            rhsCT = obj.computeElementalGradient(fGT);
            rhsS = rhsC + rhsCT;
       end        

       function rhsC = computeElementalGradient(obj,fG)
             dV    = obj.mesh.computeDvolume(obj.quadrature);
            dNdx = obj.testF.computeCartesianDerivatives(obj.quadrature);
            nDim  = size(dNdx,1);
            nNode = size(dNdx,2);
            nElem = size(dNdx,3);
            nGaus = size(dNdx,4);
            int = zeros(nNode,nElem,nDim);
            for iDime = 1:nDim
                for jDime = 1:nDim
                   % for kDime = 1:nDim
                        for igaus = 1:nGaus
                            for inode = 1:nNode
                                 fI    = squeeze(fG(iDime,:,:,igaus));
                          %      fdV    = fIJ.*dV(igaus,:);
                                 dN = squeeze(dNdx(:,inode,:,igaus))';
                                 fN    = sum(fI.*dN');
                                 intI(1,:,1) = (fN.*dV(igaus,:))';

                    %            dNK = squeeze(dNdx(kdime,inode,:,igaus))';                                
                             %   intI = dNJ.*fdV;
                                int(inode,:,iDime) = int(inode,:,iDime) + intI;
                            end
                        end
                   % end
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