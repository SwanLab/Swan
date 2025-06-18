classdef LHSIntegratorDobleWeakDivergence < handle

    properties (Access = private)
        trial
        test
        mesh
        quadrature
    end

    methods (Access = public)

        function obj = LHSIntegratorDobleWeakDivergence(cParams)
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
            dNdxTe = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTr = obj.trial.evaluateCartesianDerivatives(xV);
            dvol = obj.mesh.computeDvolume(obj.quadrature)';  % [nElem x nGauss]

            nElem  = obj.mesh.nelem;
            nDim   = size(dNdxTe,1);
            nNodeT = size(dNdxTe,2);
            nNodeR = size(dNdxTr,2);
            nGauss = size(dNdxTe,3);

            nDofTest  = nNodeT * nDim;
            nDofTrial = nNodeR * nDim;

            lhs = zeros(nDofTest, nDofTrial, nElem);

            for igaus = 1:nGauss
                for inodeT = 1:nNodeT
                    div_v = zeros(nElem,1);
                    for idim = 1:nDim
                        div_v = div_v + squeeze(dNdxTe(idim,inodeT,igaus,:));
                    end
                    for inodeR = 1:nNodeR
                        div_w = zeros(nElem,1);
                        for jdim = 1:nDim
                            div_w = div_w + squeeze(dNdxTr(jdim,inodeR,igaus,:));
                        end
                        val = div_v .* div_w .* dvol(:,igaus);
                        % Assign to all combinations of dof_test and dof_trial
                        for idim = 1:nDim
                            dof_test = (inodeT - 1)*nDim + idim;
                            for jdim = 1:nDim
                                dof_trial = (inodeR - 1)*nDim + jdim;
                                lhs(dof_test, dof_trial, :) = squeeze(lhs(dof_test, dof_trial, :)) + val;
                            end
                        end
                    end
                end
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh  = cParams.mesh;
            obj.trial = cParams.trial;
            obj.test  = cParams.test;
        end

        function createQuadrature(obj)
            orderTr = obj.trial.getOrderNum();
            orderTe = obj.test.getOrderNum();
            order = orderTr + orderTe;
            obj.quadrature = Quadrature.create(obj.mesh, order);
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun = []; % Placeholder
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
        end

    end

end
