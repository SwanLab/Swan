classdef LHSintegratorFunctionAdvection < LHSintegrator

    properties (Access = private)
        fun
    end

    methods (Access = public)

        function obj = LHSintegratorFunctionAdvection(cParams)
            obj@LHSintegrator(cParams)
            obj.fun = cParams.function;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            xV = obj.quadrature.posgp;
            shapesTest  = obj.test.computeShapeFunctions(xV);
            dNdxTr      = obj.trial.evaluateCartesianDerivatives(xV);
            fG          = squeeze(obj.fun.evaluate(xV));

            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);

            nNodeTest  = size(shapesTest,1);
            nNodeTrial = size(dNdxTr,2);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;


            lhs = zeros(nDofTest,nDofTrial,nElem);
            %             for iGaus = 1:nGaus
            %                 dV(1,1,:) = dVolu(iGaus,:)';
            %                 for iDof = 1:nDofETs
            %                     for jDof = 1:nDofETr
            %                        Ni  = shapesTest(iDof,iGaus,:);
            %                        dNj = squeeze(dNdxTr(:,jDof,:,iGaus));
            %                        df  = squeeze(fG(:,iGaus,:));
            %                        int(1,1,:) = sum(Ni*df.*dNj,1);
            %                        lhs(iDof,jDof,:) = lhs(iDof,jDof,:) + int.*dV;
            %                     end
            %
            %                 end
            %             end


            for iGaus = 1 :nGaus
                for iNode = 1:nNodeTest
                    for jNode = 1:nNodeTrial
                        for iDim = 1:obj.test.ndimf
                            for jDim = 1:obj.trial.ndimf
                                fdv = fG(iGaus,:).*dVolu(iGaus,:);
                                idof = obj.test.ndimf*(iNode-1)+iDim;
                                jdof = obj.trial.ndimf*(jNode-1)+jDim;
                                Ni = shapesTest(iNode,iGaus,:);
                                dNj = squeeze(dNdxTr(iDim,jNode,iGaus,:));
                                v = squeeze(Ni.*dNj.*(fdv'));
                                lhs(idof, jdof, :)= squeeze(lhs(idof,jdof,:)) ...
                                    + v(:);
                            end
                        end
                    end
                end
            end
        end

    end

 
end