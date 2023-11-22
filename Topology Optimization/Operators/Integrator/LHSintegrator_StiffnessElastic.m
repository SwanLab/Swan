classdef LHSintegrator_StiffnessElastic < LHSintegrator

    properties (Access = private)
        material
    end

    methods (Access = public)

        function obj = LHSintegrator_StiffnessElastic(cParams)
            obj@LHSintegrator(cParams)
            obj.material = cParams.material;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            dNdx  = obj.fun.computeCartesianDerivatives(obj.quadrature);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(obj.material.C,3);
            nNodE = size(dNdx,2);
            nDofE = nNodE*obj.fun.ndimf;
            lhs = zeros(nDofE,nDofE,nElem);
            Bcomp = obj.createBComputer(dNdx);
            Cmat = obj.material.C(:,:,:,1);
            for igaus = 1:nGaus
                Bmat = Bcomp.compute(igaus);
                dV(1,1,:) = dVolu(igaus,:)';
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
%             LHS = assembler.assemble(lhs);
            LHS = assembler.assembleFunctions(lhs, obj.fun, obj.fun);
        end
    end

end