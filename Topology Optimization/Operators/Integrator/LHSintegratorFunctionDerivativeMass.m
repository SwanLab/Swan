classdef LHSintegratorFunctionDerivativeMass < handle

    properties (Access = private)
        mesh
        test, trial
        quadrature
        quadratureOrder
        fun
    end

    methods (Access = public)

        function obj = LHSintegratorFunctionDerivativeMass(cParams)
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
            quad = obj.quadrature;
            xV = quad.posgp;
            shapesTest  = obj.test.computeShapeFunctions(xV);
            shapesTrial = obj.trial.computeShapeFunctions(xV);
            dVolu  = obj.mesh.computeDvolume(quad);
            nGaus  = quad.ngaus;
            nElem  = size(dVolu,2);

            nNodeTest  = size(shapesTest,1);
            nNodeTrial = size(shapesTrial,1);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;

            dfdx       = obj.fun.evaluateGradient(xV);

            M = zeros(nDofTest, nDofTrial, nElem);
            % lhs = zeros(nDofTest/2, nDofTrial/2, nElem);
            % for igaus = 1:nGaus
            %     dv(1,1,:) = dVolu(igaus,:);
            %     Nv = shapesTest(:,igaus);
            %     Nu = shapesTrial(:,igaus);
            %     NvNu = Nv*Nu';
            %     Aij = bsxfun(@times,NvNu,dv);
            %     lhs = lhs + Aij;
            % end
            for iGaus = 1 :nGaus
                for iNode = 1:nNodeTest
                    for jNode = 1:nNodeTrial
                        for iDim = 1:obj.test.ndimf
                            for jDim = 1:obj.trial.ndimf
                                fGI = squeezeParticular(dfdx.fValues(iDim,iGaus,:),1);
                                fdv = fGI.*dVolu(iGaus,:);
                                idof = obj.test.ndimf*(iNode-1)+iDim;
                                jdof = obj.trial.ndimf*(jNode-1)+jDim;
                                Ni = shapesTest(iNode,iGaus,:);
                                Nj = shapesTrial(jNode,iGaus,:);
                                v = squeeze(Ni.*Nj.*fdv);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:);
                            end
                        end
                    end
                end
            end
            lhs = M;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.fun = cParams.function;
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