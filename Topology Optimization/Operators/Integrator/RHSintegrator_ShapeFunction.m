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
        end

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
            fG   = fun.evaluate(xV);
            dV   = obj.mesh.computeDvolume(quad);
            N = test.computeShapeFunctions(xV);
            nDofElem  = size(N,1);
            nElem     = obj.mesh.nelem;
            nGaus     = quad.ngaus;
            nFlds     = size(fG,1);
            int = zeros(nDofElem,nElem);
            for iField = 1:nFlds
                for iDof = 1:nDofElem
                    for iGaus = 1:nGaus
                        dVg(:,1) = dV(iGaus, :);
                        fV   = squeeze(fG(iField,iGaus,:));
                        Ni   = squeeze(N(iDof,iGaus,:));
                        fNdV(1,:) = Ni.*fV.*dVg;
                        int(iDof,:) = int(iDof,:) + fNdV;
                    end
                end
            end
            rhsC = pagetranspose(int);
        end

        function f = assembleIntegrand(obj,test,rhsElem)
            integrand = rhsElem;
            connec   = test.computeDofConnectivity()';
            ndofs    = max(max(connec));
            nDofElem = size(connec,2);
            f = zeros(ndofs,1);
            for iDof = 1:nDofElem
                int = integrand(:,iDof);
                con = connec(:,iDof);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadType);
            obj.quadrature = q;
        end

    end

end
