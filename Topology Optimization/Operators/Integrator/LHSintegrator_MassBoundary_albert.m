classdef LHSintegrator_MassBoundary_albert < handle

    properties (Access = private)
        boundaryMeshJoined
        localGlobalConnecBd
        quadType
        quadrature
        nnodes
    end

    methods (Access = public)

        function obj = LHSintegrator_MassBoundary_albert(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj,trial,test)
            obj.createQuadrature(trial,test)
            LHS = obj.computeBoundaryMassMatrix(trial,test);
        end

    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.quadType = cParams.quadType;
            obj.boundaryMeshJoined  = cParams.boundaryMeshJoined;
            obj.localGlobalConnecBd = cParams.localGlobalConnecBd;
            obj.nnodes              = cParams.nnodes;
        end

        function Mr = computeBoundaryMassMatrix(obj,trial,test)
            nDofs = obj.nnodes*trial.ndimf;
            lhsg = sparse(nDofs,trial.nDofs);
            a.type = 'MassMatrix';
            a.mesh = obj.boundaryMeshJoined;

            a.test  = test;
            a.trial = trial;


            lhs = LHSIntegrator.create(a);

            lhs = lhs.compute();
            [iLoc,jLoc,vals] = find(lhs);
            l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
            l2g_dof = l2g_dof(:);
            jGlob = l2g_dof(jLoc);
            iGlob = l2g_dof(iLoc);
            lhsg = lhsg + sparse(iGlob,jLoc,vals, nDofs, trial.nDofs);

            % local2global(sL.mesh.connec(:)) = sL.bMesh.globalConnec(:);
            % [iLoc,jLoc,vals] = find(LHS); % !!! iLoc, jLoc should come from P1Fun
            % iGlob = local2global(iLoc);
            % jGlob = local2global(jLoc);

            %         LHSadd = sparse(iGlob,jGlob,vals, ndof, ndof);
            %         LHSg = LHSg + LHSadd;
            Mr = lhsg;
        end

        function createQuadrature(obj,fun,test)
            if isempty(obj.quadType)
                orderTr = fun.getOrderNum();
                orderTe = test.getOrderNum();
                order = orderTr + orderTe;
            else
                order = obj.quadType;
            end
            q = Quadrature.create(obj.boundaryMeshJoined,order);
            obj.quadrature = q;
        end


    end

end