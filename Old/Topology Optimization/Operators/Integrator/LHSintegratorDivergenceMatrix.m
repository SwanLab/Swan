classdef LHSintegratorDivergenceMatrix < handle

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSintegratorDivergenceMatrix(cParams)
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
            dNdxTe  = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTr  = obj.trial.evaluateCartesianDerivatives(xV);

            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);

            nNodeTest  = size(dNdxTe,2);
            nNodeTrial = size(dNdxTr,2);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;


            lhs = zeros(nDofTest,nDofTrial,nElem);

            for iGaus = 1 :nGaus
                for iNode = 1:nNodeTest
                    for jNode = 1:nNodeTrial
                        for iDim = 1:obj.test.ndimf
                            for jDim = 1:obj.trial.ndimf
                                dVG   = dVolu(iGaus,:)';
                                idof = obj.test.ndimf*(iNode-1)+iDim;
                                jdof = obj.trial.ndimf*(jNode-1)+jDim;
                                dNi = squeeze(dNdxTr(iDim,iNode,:,iGaus));
                                dNj = squeeze(dNdxTr(jDim,jNode,:,iGaus));
                                v = squeeze(dNi.*dNj.*dVG);
                                lhs(idof, jdof, :)= squeeze(lhs(idof,jdof,:)) ...
                                    + v(:);
                            end
                        end
                    end
                end
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

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
        end

    end

end