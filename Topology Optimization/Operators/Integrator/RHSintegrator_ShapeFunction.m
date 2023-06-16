classdef RHSintegrator_ShapeFunction < RHSintegrator

    properties (Access = private)
        test
        quadratureOrder
    end

    methods (Access = public)

        function obj = RHSintegrator_ShapeFunction(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function rhs = compute(obj, fun)
            rhsElem = obj.computeElementalRHS(fun);
            rhs = obj.assembleIntegrand(fun,rhsElem);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.test = cParams.test;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = obj.test.order;
            end
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end        

        function rhsC = computeElementalRHS(obj, fun)
            quad = obj.quadrature;
            fG     = fun.evaluate(quad.posgp); 
            dV     = obj.computeDvolume();
            N = obj.test.computeShapeFunctions(quad);
            N = permute(N, [1 3 2]);
            nNode  = size(N,1);
            nElem  = obj.mesh.nelem;
            nGaus  = quad.ngaus;
            shapes = repmat(N, [1 1 nElem]);
            int = zeros(nNode,nElem,fun.ndimf);
            for iDim = 1:fun.ndimf
                fGI = squeeze(fG(iDim,:,:));
                for iGaus = 1:nGaus
                    fdv = fGI(iGaus,:).*dV(iGaus,:);
                    shape = shapes(:, :, iGaus);
                    int = int + bsxfun(@times,shape,fdv);
                end
            end
            rhsC = permute(int,[2 1 3]);
        end

        function f = assembleIntegrand(obj,fun,rhsElem)
            integrand = rhsElem;
            connec = obj.test.computeDofConnectivity()';
            ndofs = max(max(connec));
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for iNode = 1:nnode
                int = integrand(:,iNode);
                con = connec(:,iNode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end

        function dV = computeDvolume(obj)
            q = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);
        end

    end

end
