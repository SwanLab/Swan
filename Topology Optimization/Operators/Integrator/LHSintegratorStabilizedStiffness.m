classdef LHSintegratorStabilizedStiffness < handle %LHSintegrator

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
        Tau
    end

    methods (Access = public)

        function obj = LHSintegratorStabilizedStiffness(cParams)
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

            t   = obj.Tau;
            lhs = zeros(nDofETs,nDofETr,nElem);
            for igaus = 1:nGaus
                BmatTs = squeezeParticular(dNdxTs(:,:,igaus,:),3);
                BmatTr = squeezeParticular(dNdxTr(:,:,igaus,:),3);
                dV(1,1,:) = dVolu(igaus,:)';
                Bt   = permute(BmatTs,[2 1 3]);
                BtC  = pagemtimes(t,Bt);
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
            obj.Tau   = cParams.Tau;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = obj.trial.orderTextual;
            end
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun     = [];
            assembler = AssemblerFun(s);
            LHS       = assembler.assemble(lhs, obj.test, obj.trial);
        end

    end

end