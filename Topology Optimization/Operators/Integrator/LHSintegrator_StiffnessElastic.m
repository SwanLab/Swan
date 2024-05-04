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
            xV = obj.quadrature.posgp;
            dNdx  = obj.fun.evaluateCartesianDerivatives(xV);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = obj.mesh.nelem;
            nNodE = size(dNdx,2);
            nDofE = nNodE*obj.fun.ndimf;
            lhs = zeros(nDofE,nDofE,nElem);
            Bcomp = obj.createBComputer(dNdx);
            Cmat  = obj.material.evaluate(xV);
            %Cmat = obj.material;
            for igaus = 1:nGaus
                Bmat = Bcomp.compute(igaus);
                C    = squeeze(Cmat(:,:,igaus,:));
                dV(1,1,:) = dVolu(igaus,:)';
                Bt   = permute(Bmat,[2 1 3]);
                BtC  = pagemtimes(Bt,C);
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
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.fun, obj.fun);
        end
    end

end