classdef LHSintegrator_MassDomain < handle

    properties (Access = private)
        mesh
        test, trial
        domain
        quadrature
        quadratureOrder
        l2g, ndofs
    end

    methods (Access = public)

        function obj = LHSintegrator_MassDomain(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.setQuadratureOrder(cParams);

            obj.ndofs = numel(obj.trial.fValues);
            if isfield(cParams, 'domain')
                obj.domain = cParams.domain;
                [obj.trial, obj.mesh, obj.l2g] = obj.trial.evaluateBoundarySubdomain(obj.domain);
                obj.test = obj.test.evaluateBoundarySubdomain(obj.domain);
            else
                obj.domain = [];
            end
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

        function lhs = computeElementalLHS(obj)
            % obj.ndofs = numel(obj.trial.fValues);
            % if ~isempty(obj.domain)
            %     [obj.trial, obj.mesh, obj.l2g] = obj.trial.evaluateBoundarySubdomain(obj.domain);
            %     obj.test = obj.test.evaluateBoundarySubdomain(obj.domain);
            % end

            quad = obj.quadrature;
            shapesTest  = obj.test.computeShapeFunctions(quad);
            shapesTrial = obj.trial.computeShapeFunctions(quad);
            dVolu  = obj.mesh.computeDvolume(quad);

            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);
            nNodeTest  = size(shapesTest,1);
            nNodeTrial = size(shapesTrial,1);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;

            M = zeros(nDofTest, nDofTrial, nElem);
%             for igaus = 1:nGaus
%                 dv(1,1,:) = dVolu(igaus,:);
%                 Nv = shapesTest(:,igaus);
%                 Nu = shapesTrial(:,igaus);
%                 NvNu = Nv*Nu';
%                 Aij = bsxfun(@times,NvNu,dv);
%                 lhs = lhs + Aij;
%             end
            for igauss = 1 :nGaus
                for inode= 1:nNodeTest
                    for jnode= 1:nNodeTrial
                        for iunkn= 1:obj.test.ndimf
                       %     for junkn= 1:obj.trial.ndimf
                                idof = obj.test.ndimf*(inode-1)+iunkn;
                                jdof = obj.trial.ndimf*(jnode-1)+iunkn;
                                dvol = dVolu(igauss,:);
                                Ni = shapesTest(inode,igauss,:);
                                Nj = shapesTrial(jnode,igauss,:);
                                v = squeeze(Ni.*Nj);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:).*dvol';
                       %     end
                        end
                    end
                end
            end
            lhs = M;

        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
            [iLoc,jLoc,vals] = find(LHS);
            globalDofConnec = ((obj.l2g*obj.test.ndimf)' - (1:-1:0))';
            l2g_dof = globalDofConnec(:);
            iGlob = l2g_dof(iLoc);
            jGlob = l2g_dof(jLoc);
            nDofs = obj.ndofs;
            LHS = sparse(iGlob,jGlob,vals, nDofs, nDofs);
        end

    end

end