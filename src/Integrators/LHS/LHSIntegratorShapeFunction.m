classdef LHSIntegratorShapeFunction < handle

    properties (Access = private)
        quadType
        mesh
        quadrature
    end

    methods (Access = public)
        function obj = LHSIntegratorShapeFunction(cParams)
            obj.init(cParams);
        end


        function lhs = compute(obj,fun,test)
            obj.createQuadrature(fun,test);
            lhsElem = obj.computeElementalLHS(fun,test);
            lhs = obj.assembleIntegrand(test,fun,lhsElem);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadType = cParams.quadType;
            obj.mesh     = cParams.mesh;
        end

        function int = computeElementalLHS(obj,fun,test)
            quad = obj.quadrature;
            xV   = quad.posgp;
            fG   = fun.evaluate(xV);
%             fG   = squeezeParticular(fG,2);
            dV   = obj.mesh.computeDvolume(quad);
            N = test.computeShapeFunctions(xV);
            nNodeElem  = size(N,1);
            nElem     = obj.mesh.nelem;
            nGaus     = quad.ngaus;
            nFlds     = size(fG,1);
            nDofElem = nNodeElem*nFlds;
            int = zeros(nDofElem,nFlds,nElem);
            for iField = 1:nFlds
                for iNode = 1:nNodeElem
                    for iGaus = 1:nGaus
                        dVg(:,1) = dV(iGaus, :);
                        fV   = squeeze(fG(iField,iGaus,:));
                        Ni   = squeeze(N(iNode,iGaus,:));
                        fNdV(1,:) = Ni.*fV.*dVg;
                        iDof = nFlds*(iNode-1) + iField;
                       int(iDof,iField,:) = squeezeParticular(int(iDof,iField,:),1) + fNdV;
                    end
                end
            end
            
        end

        function C = assembleIntegrand(obj,test,fun,lhsElem)
            integrand = lhsElem;
            nFields = fun.ndimf;
            connec   = test.getDofConnec();
            ndofs    = max(max(connec));
            nDofElem = size(connec,2);
            C = zeros(ndofs,nFields);
            for iField = 1:nFields
                for iDof = 1:nDofElem
                    int = squeezeParticular(integrand(iDof,iField,:),1);
                    con = connec(:,iDof);
                    C(:,iField) = C(:,iField) + accumarray(con,int,[ndofs,1],@sum,0);
                end
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
