classdef LHSintegrator_Stiffness < handle %LHSintegrator

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSintegrator_Stiffness(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
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

            BcompTs = obj.createBComputer(obj.test, dNdxTs);
            BcompTr = obj.createBComputer(obj.trial, dNdxTr);
            lhs = zeros(nDofETs,nDofETr,nElem);
            for igaus = 1:nGaus
                BmatTs = BcompTs.compute(igaus);
                BmatTr = BcompTr.compute(igaus);
                dV(1,1,:) = dVolu(igaus,:)';
                Bt   = permute(BmatTs,[2 1 3]);
                BtCB = pagemtimes(Bt, BmatTr);
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
            LHS = assembler.assembleFunctions(lhs, obj.test, obj.trial);
        end

    end

end