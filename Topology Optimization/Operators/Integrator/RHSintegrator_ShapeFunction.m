classdef RHSintegrator_ShapeFunction < handle

    properties (Access = private)
        quadType
        mesh
        quadrature
    end

    methods (Access = public)
        function obj = RHSintegrator_ShapeFunction(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end


        function rhs = compute(obj,fun,test)
            rhsElem = obj.computeElementalRHS(fun,test);
            rhs = obj.assembleIntegrand(test,rhsElem);

%        function RHS = compute(obj,fun,test)
%            quad = obj.quadrature;
%            xV   = quad.posgp;
%            dV   = obj.mesh.computeDvolume(quad);
%            shapes = test.computeShapeFunctions(quad);
%            nGaus = quad.ngaus;
%            nFlds = fun.ndimf;
%            nDofElem = size(shapes,1);
%            conne = test.computeDofConnectivity';
%            nDofs = max(conne,[],"all");
%            fGaus = fun.evaluate(xV);
%            f     = zeros(nDofs,nFlds);
%            for iField = 1:nFlds
%                for igaus = 1:nGaus
%                    dVg(:,1) = dV(igaus, :);
%                    fG = squeeze(fGaus(iField,igaus,:));
%                    for idof = 1:nDofElem
%                        dofs = conne(:,idof);
%                        Ni = shapes(idof,igaus);
%                        int = Ni*fG.*dVg;
%                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
%                    end
%                end
%            end
%            RHS = f;
%
%        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadType = cParams.quadType;
            obj.mesh     = cParams.mesh;
        end

        function rhsC = computeElementalRHS(obj,fun,test)
            quad = obj.quadrature;
            xV   = quad.posgp;
            fG     = obj.fun.evaluate(xV);
            dV     = obj.computeDvolume();
            N = test.computeShapeFunctions(quad);
            N = permute(N, [1 3 2]);
            nNode  = size(N,1);
            nElem  = obj.mesh.nelem;
            nGaus  = quad.ngaus;
	    nFlds  = fun.ndimf;
            shapes = repmat(N, [1 1 nElem]);
            int = zeros(nNode,nElem);
            for iField = 1:nFlds
            	for iGaus = 1:nGaus
            	    fV  = squeeze(fG(iField,igaus,:));
            	    fdv = fV.*dV(iGaus,:);
               	    shape = shapes(:, :, iGaus);
                    int = int + bsxfun(@times,shape,fdv);
	        end
            end
            rhsC = transpose(int);
        end

        function f = assembleIntegrand(obj,test,rhsElem)
            integrand = rhsElem;
            connec = test.computeDofConnectivity()';
            ndofs  = max(max(connec));
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

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadType);
            obj.quadrature = q;

        end

    end

end
