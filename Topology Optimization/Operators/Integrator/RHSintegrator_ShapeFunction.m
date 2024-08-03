classdef RHSintegrator_ShapeFunction < handle

    properties (Access = private)
        quadType
        mesh
        quadrature
    end

    methods (Access = public)
        function obj = RHSintegrator_ShapeFunction(cParams)
            obj.init(cParams);
        end


        function rhs = compute(obj,fun,test)
            obj.createQuadrature(fun,test);
            rhsElem = obj.computeElementalRHS(fun,test);
            rhs = obj.assembleIntegrand(test,rhsElem);
        end

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
            nNodeElem  = size(N,1);
            nElem     = obj.mesh.nelem;
            nGaus     = quad.ngaus;
            nFlds     = size(fG,1);
            nDofElem = nNodeElem*nFlds;
            int = zeros(nDofElem,nElem);
            for iField = 1:nFlds
                for iNode = 1:nNodeElem
                    for iGaus = 1:nGaus
                        dVg(:,1) = dV(iGaus, :);
                        fV   = squeeze(fG(iField,iGaus,:));
                        Ni   = squeeze(N(iNode,iGaus,:));
                        fNdV(1,:) = Ni.*fV.*dVg;
                        iDof = nFlds*(iNode-1) + iField;
                        int(iDof,:) = int(iDof,:) + fNdV;
                    end
                end
            end
            rhsC = transpose(int);
        end

        function f = assembleIntegrand(obj,test,rhsElem)
            integrand = rhsElem;
            connec   = test.getConnec();
            ndofs    = max(max(connec));
            nDofElem = size(connec,2);
            f = zeros(ndofs,1);
            for iDof = 1:nDofElem
                int = integrand(:,iDof);
                con = connec(:,iDof);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end

        function createQuadrature(obj,fun,test)
            if isempty(obj.quadType)
                orderTr = fun.getOrderNum();
                orderTe = test.getOrderNum();
                order = orderTr + orderTe;
            else
                order = obj.quadType;
            end
            q = Quadrature.create(obj.mesh,order);
            obj.quadrature = q;
        end

    end

end
