classdef LHSintegrator_Stiffness < LHSintegrator

    methods (Access = public)

        function obj = LHSintegrator_Stiffness(cParams)
            obj@LHSintegrator(cParams)
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            xV = obj.quadrature.posgp;
            dNdxTs = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTr = obj.trial.evaluateCartesianDerivatives(xV);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);

            nNodETs = size(dNdxTs,2);
            nDofETs = nNodETs*obj.test.ndimf;
            nNodETr = size(dNdxTr,2);
            nDofETr = nNodETr*obj.trial.ndimf;

            lhs = zeros(nDofETs,nDofETr,nElem);
            for igaus = 1:nGaus
                BmatTs = squeezeParticular(dNdxTs(:,:,igaus,:),3);
                BmatTr = squeezeParticular(dNdxTr(:,:,igaus,:),3);
                dV(1,1,:) = dVolu(igaus,:)';
                Bt   = permute(BmatTs,[2 1 3]);
                BtCB = pagemtimes(Bt, BmatTr);
                lhs = lhs + bsxfun(@times, BtCB, dV);
            end
        end

    end

    methods (Access = private)

        function Bcomp = createBComputer(obj, fun, dNdx)
            s.fun  = fun;
            s.dNdx = dNdx;
            Bcomp = BMatrixComputer(s);
        end

    end

end