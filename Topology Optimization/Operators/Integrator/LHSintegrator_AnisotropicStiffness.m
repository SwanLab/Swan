classdef LHSintegrator_AnisotropicStiffness < handle %LHSintegrator

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder

        CAnisotropic
        Celas
        alphaDeg
    end

    methods (Access = public)

        function obj = LHSintegrator_AnisotropicStiffness(cParams)
            obj.init(cParams);
            obj.createQuadrature();
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
            dNdxTs = obj.test.computeCartesianDerivatives(obj.quadrature);
            dNdxTr = obj.trial.computeCartesianDerivatives(obj.quadrature);
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

        function init(obj, cParams)
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = obj.trial.order;
            end
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function Bcomp = createBComputer(obj, fun, dNdx)
            s.fun  = fun;
            s.dNdx = dNdx;
            Bcomp = BMatrixComputer(s);
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
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