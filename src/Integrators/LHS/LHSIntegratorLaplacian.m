classdef LHSIntegratorLaplacian < handle

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSIntegratorLaplacian(cParams)
            obj.initLaplacian(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function LHSe = computeElementalLHS(obj)
            xV = obj.quadrature.posgp;
            shapesTs = obj.test.evaluateCartesianDerivatives(xV);
            shapesTr = obj.trial.evaluateCartesianDerivatives(xV);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);

            nGaus = obj.quadrature.ngaus;
            nNodETs = size(shapesTs,2);
            nDofETs = nNodETs*obj.test.ndimf;
            nNodETr = size(shapesTr,2);
            nDofETr = nNodETr*obj.trial.ndimf;
            nElem = size(dVolu,2);
            
            lhs = zeros(nDofETs, nDofETr, nElem);
            for iGaus = 1:nGaus
                dNdxTs = squeeze(shapesTs(:,:,iGaus,:));
                dNdxTr = squeeze(shapesTr(:,:,iGaus,:));
                dNTs = obj.reshapeGradient(dNdxTs);
                dNTr = obj.reshapeGradient(dNdxTr);
                dV(1,1,:) = dVolu(iGaus,:)';
                Bt   = permute(dNTs,[2 1 3]);
                BtB = pagemtimes(Bt,dNTr);
                lhs = lhs + pagemtimes(BtB,dV);
            end
            LHSe = lhs;
        end

    end

    methods (Access = private)

        function initLaplacian(obj, cParams)
            obj.mesh  = cParams.mesh;
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
        end

        function createQuadrature(obj)
            orderTr = obj.trial.getOrderNum();
            orderTe = obj.test.getOrderNum();
            order = orderTr + orderTe;
            quad = Quadrature.create(obj.mesh,order);
            obj.quadrature = quad;
        end

        function dN = reshapeGradient(obj,dNdx)
            nunkn = size(dNdx,1);
            nnode = size(dNdx,2);
            nelem = size(dNdx,3);
            dN = zeros(4,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                dN(1,j,:)   = dNdx(1,i,:);
                dN(2,j+1,:) = dNdx(1,i,:);
                dN(3,j,:)   = dNdx(2,i,:);
                dN(4,j+1,:) = dNdx(2,i,:);
            end
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
        end

    end

end