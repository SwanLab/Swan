classdef LHSIntegratorAnisotropicStiffness < LHSIntegrator

    properties (Access = private)
        CAnisotropic
        Celas
        alphaDeg
    end

    methods (Access = public)

        function obj = LHSIntegratorAnisotropicStiffness(cParams)
            obj@LHSIntegrator(cParams)
            obj.initAnisotropicTensor(cParams);
        end

        function LHS = compute(obj)
            obj.assemblyCMatrix();
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

            Cmat    = obj.Celas;
            BcompTs = obj.createBComputer(obj.test, dNdxTs);
            BcompTr = obj.createBComputer(obj.trial, dNdxTr);
            lhs = zeros(nDofETs,nDofETr,nElem);
            for iGaus = 1:nGaus
                BmatTs = BcompTs.compute(iGaus);
                BmatTr = BcompTr.compute(iGaus);
                dV(1,1,:) = dVolu(iGaus,:)';
                Bt   = permute(BmatTs,[2 1 3]);
                BtC  = pagemtimes(Bt,Cmat);
                BtCB = pagemtimes(BtC, BmatTr);
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

        function initAnisotropicTensor(obj,cParams)
            CLocal = cParams.CAnisotropic;
            obj.alphaDeg = cParams.aniAlphaDeg;
            obj.CAnisotropic = obj.rotateAnisotropicMatrix(CLocal);
        end

        function CGlobal = rotateAnisotropicMatrix(obj,CLocal)
            R = [cosd(obj.alphaDeg),-sind(obj.alphaDeg)
                sind(obj.alphaDeg), cosd(obj.alphaDeg)];
            CGlobal = R*CLocal*R';
        end

        function assemblyCMatrix(obj)
            nelem = size(obj.mesh.connec,1);
            C = repmat(obj.CAnisotropic, [1 1 nelem]);
            obj.Celas = C;
        end
    end

end