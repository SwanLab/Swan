classdef LHSintegratorAnisotropicStiffness < LHSintegrator

    properties (Access = private)
        CAnisotropic
        Celas
        alphaDeg
    end

    methods (Access = public)

        function obj = LHSintegratorAnisotropicStiffness(cParams)
            obj@LHSintegrator(cParams);
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
            dNdx  = obj.fun.computeCartesianDerivatives(obj.quadrature);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);
            nNodE = size(dNdx,2);
            nDofE = nNodE*obj.fun.ndimf;
            lhs = zeros(nDofE,nDofE,nElem);
            Cmat   = obj.Celas;
            Bcomp = obj.createBComputer(dNdx);
            for iGaus = 1:nGaus
                Bmat = Bcomp.compute(iGaus);
                dV(1,1,:) = dVolu(iGaus,:)';
                Bt   = permute(Bmat,[2 1 3]);
                BtC  = pagemtimes(Bt,Cmat);
                BtCB = pagemtimes(BtC, Bmat);
                lhs = lhs + bsxfun(@times, BtCB, dV);
            end
        end

    end

    methods (Access = private)

        function Bcomp = createBComputer(obj, dNdx)
            s.fun  = obj.fun;
            s.dNdx = dNdx;
            Bcomp = BMatrixComputer(s);
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = obj.fun; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs);
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