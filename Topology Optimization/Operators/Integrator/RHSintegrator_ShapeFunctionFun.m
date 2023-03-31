classdef RHSintegrator_ShapeFunctionFun < handle

    properties (Access = private)
        mesh
        quadrature
    end

    methods (Access = public)

        function obj = RHSintegrator_ShapeFunctionFun(cParams)
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
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        function rhsC = computeElementalRHS(obj, fun)
            quad = obj.quadrature;
%             fG     = obj.fun.evaluate(quad.posgp);
            fG     = squeeze(fun.evaluate(quad.posgp)); % only used in 1d funs so far
            dV     = obj.computeDvolume();
            N = fun.computeShapeFunctions(quad);
            N = permute(N, [1 3 2]);
            nNode  = size(N,1);
            nElem  = obj.mesh.nelem;
            nGaus  = quad.ngaus;
            shapes = repmat(N, [1 1 nElem]);
            int = zeros(nNode,nElem);
            for iGaus = 1:nGaus
                fdv = fG(iGaus,:).*dV(iGaus,:);
                shape = shapes(:, :, iGaus);
                int = int + bsxfun(@times,shape,fdv);
            end
            rhsC = transpose(int);
        end

        function f = assembleIntegrand(obj,fun,rhsElem)
            integrand = rhsElem;
            connec = fun.computeDofConnectivity()';
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
