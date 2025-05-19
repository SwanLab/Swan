classdef LHSIntegratorLaplacian < handle

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
        material
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

%             material = obj.material;
            
            Cmat = obj.material.mu;
            
            lhs = zeros(nDofETs, nDofETr, nElem);
            for iGaus = 1:nGaus
                dNdxTs = squeeze(shapesTs(:,:,iGaus,:));
                dNdxTr = squeeze(shapesTr(:,:,iGaus,:));
                BmatTs = obj.computeB(dNdxTs);
                BmatTr = obj.computeB(dNdxTr);
                dV(1,1,:) = dVolu(iGaus,:)';
                Bt   = permute(BmatTs,[2 1 3]);
                BtC  = pagemtimes(Bt,Cmat);
                BtCB = pagemtimes(BtC, BmatTr);
                %BtB = pagemtimes(Bt,BmatTr);
                lhs = lhs + bsxfun(@times, BtCB, dV);
            end
            LHSe = lhs;
        end

    end

    methods (Access = private)

        function initLaplacian(obj, cParams)
            obj.mesh  = cParams.mesh;
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.material = cParams.material;
        end

        function createQuadrature(obj)
            orderTr = obj.trial.getOrderNum();
            orderTe = obj.test.getOrderNum();
            order = orderTr + orderTe;
            quad = Quadrature.create(obj.mesh,order);
            obj.quadrature = quad;
        end

        function B = computeB(obj,dNdx)
            nunkn = size(dNdx,1);
            nnode = size(dNdx,2);
            nelem = size(dNdx,3);
            B = zeros(4,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)   = dNdx(1,i,:);
                B(2,j+1,:) = dNdx(1,i,:);
                B(3,j,:)   = dNdx(2,i,:);
                B(4,j+1,:) = dNdx(2,i,:);
            end
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
        end

    end

end