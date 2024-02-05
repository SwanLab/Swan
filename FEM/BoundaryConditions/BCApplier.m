classdef BCApplier < handle
    
    % Goal: group dirichlet and neumann conditions
    % to allow multiple boundary conditions at the same time
    % Use: BCApplier.computeLinearConditionsMatrix()

    properties (Access = public)
    end
    
    properties (Access = private)
        mesh

        dirichlet_dofs, dirichlet_vals
        dirichletFun
        pointload_dofs, pointload_vals
        pointloadFun
        periodic_leader, periodic_follower
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = BCApplier(cParams)
            obj.init(cParams)
        end
        
        function Ct = computeLinearConditionsMatrix(obj)
            % generalize lagrange multiplier -> dirac and stuff
            dir_dofs = obj.dirichlet_dofs;
            nDofs = obj.dirichletFun.nDofs;
            nDirich = length(dir_dofs);
            Ct = sparse(1:nDirich, dir_dofs, 1, nDirich, nDofs);
        end

        function Ct = computeP1LinearConditionsMatrix(obj)
            nDofs = obj.dirichletFun.nDofs;
            mesh_left = obj.mesh.createBoundaryMesh{1};
            % Check LHSintegrator_MassBoundary
            a.mesh = mesh_left.mesh;
            a.test = P1Function.create(mesh_left.mesh, 2); % from elastic
            a.trial = P1Function.create(mesh_left.mesh, 2); % dLambda
            a.type = 'MassMatrix';
            lhs = LHSintegrator.create(a);
            LHS = lhs.compute();

            local_dofConnec = obj.computeDofConnectivity(mesh_left.mesh.connec);
            global_dofConnec = obj.computeDofConnectivity(mesh_left.globalConnec);

            local2global(local_dofConnec(:)) = global_dofConnec(:);
            [iLoc,jLoc,vals] = find(LHS); % !!! iLoc, jLoc should come from P1Fun
            iGlob = local2global(iLoc);
            jGlob = local2global(jLoc);

            Ct = sparse(iGlob,jGlob,vals, nDofs, nDofs);
            Ct = Ct(unique(iGlob), :);

            % Check LHSintegrator_MassBoundary
            isLeft =  @(coor) (abs(coor(:,1)) == 0);
            b.mesh = obj.mesh;
            b.test = P1Function.create(obj.mesh, 2); % from elastic
            b.trial = P1Function.create(obj.mesh, 2); % dLambda
            b.domain = isLeft;
            b.type = 'MassDomain';
            lhs3 = LHSintegrator.create(b);
            LHS3 = lhs3.compute();

            [aa,~] = find(LHS3);
            Ct2 = LHS3(unique(aa),:);
        end

        function Ct = computeP0LinearConditionsMatrix(obj)
            nDofs = obj.dirichletFun.nDofs;
            mesh_left = obj.mesh.createBoundaryMesh{1};
            % Check LHSintegrator_MassBoundary
            a.mesh = mesh_left.mesh;
            a.test  = P0Function.create(mesh_left.mesh, 2); % from elastic
            a.trial = P1Function.create(mesh_left.mesh, 2); % dLambda
            a.type = 'MassMatrix';
            lhs = LHSintegrator.create(a);
            LHS = lhs.compute();

            local_dofConnec = obj.computeDofConnectivity(mesh_left.mesh.connec);
            global_dofConnec = obj.computeDofConnectivity(mesh_left.globalConnec);

            local2global(local_dofConnec(:)) = global_dofConnec(:);
            [iLoc,jLoc,vals] = find(LHS); % !!! iLoc, jLoc should come from P1Fun
            iGlob = local2global(iLoc);
            jGlob = local2global(jLoc);

            Ct = sparse(iGlob,jGlob,vals, nDofs, nDofs);
            Ct = Ct(unique(iGlob), :);

            % Check LHSintegrator_MassBoundary
            isLeft =  @(coor) (abs(coor(:,1)) == 0);
            b.mesh = obj.mesh;
            b.test  = P0Function.create(obj.mesh, 2); % from elastic
            b.trial = P1Function.create(obj.mesh, 2); % dLambda
            b.domain = isLeft;
            b.type = 'MassDomain';
            lhs3 = LHSintegrator.create(b);
            LHS3 = lhs3.compute();

            [aa,~] = find(LHS3);
            Ct = LHS3(unique(aa),:);
        end

        function dofConnec = computeDofConnectivity(obj, conne)
            nDimf  = 2;
            nNode  = size(conne, 2);
            nDofsE = nNode*nDimf;
            dofsElem  = zeros(nDofsE,size(conne,1));
            for iNode = 1:nNode
                for iUnkn = 1:nDimf
                    idofElem   = nDimf*(iNode - 1) + iUnkn;
                    globalNode = conne(:,iNode);
                    idofGlobal = nDimf*(globalNode - 1) + iUnkn;
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            dofConnec = dofsElem;
        end

        function Ct = computeLinearPeriodicConditionsMatrix(obj)
            per_lead = obj.periodic_leader;
            per_fllw = obj.periodic_follower;
            nDofs = obj.dirichletFun.nDofs;
            nPer = length(per_lead);
            Ct = sparse([(1:nPer)', (1:nPer)'], [per_lead, per_fllw], [ones(size(per_lead,1),1), -ones(size(per_lead,1),1)], nPer, nDofs); % !!
        end

        function Ct = computeSingleDirichletPeriodicCondition(obj, iVoigt, nVoigt)
            per_lead = obj.periodic_leader;
            per_fllw = obj.periodic_follower;
            nDofs = obj.dirichletFun.nDofs;
            % Dirichlet stuff
            basis   = diag(ones(nVoigt,1));
            vstrain = basis(iVoigt,:);
            s.mesh           = obj.mesh;
            s.dirichlet_dofs = obj.dirichlet_dofs;
            s.ndofs          = nDofs;
            s.vstrain        = vstrain;
            comp = MicroDispMonolithicBCApplier(s);
            [CtDir, CtDirPer] = comp.getLHSMatrix();
            % Periodic stuff
            nDofsPerBorder = length(obj.periodic_leader)/4; % 4 because 2D
            xx_bottom = 1:nDofsPerBorder;
            xx_left   = nDofsPerBorder+1 : 2*nDofsPerBorder;
            yy_bottom = 2*nDofsPerBorder + 1 : 3*nDofsPerBorder;
            yy_left   = 3*nDofsPerBorder + 1 : 4*nDofsPerBorder;
            CtPer = sparse([(1:nDofsPerBorder)', (1:nDofsPerBorder)'; ... % xx
                         (nDofsPerBorder+1:2*nDofsPerBorder)', (nDofsPerBorder+1:2*nDofsPerBorder)'; ... % xy
                         (nDofsPerBorder+1:2*nDofsPerBorder)', (nDofsPerBorder+1:2*nDofsPerBorder)'; ... % xy
                         (2*nDofsPerBorder+1:3*nDofsPerBorder)', (2*nDofsPerBorder+1:3*nDofsPerBorder)'; ... % yy
                         ], ...
                         [per_lead(xx_bottom), per_fllw(xx_bottom); ... % xx
                         per_lead(xx_left), per_fllw(xx_left); ... % xy
                         per_lead(yy_bottom), per_fllw(yy_bottom); ... % xy
                         per_lead(yy_left), per_fllw(yy_left) % yy
                         ], ...
                         [ones(length(per_lead),1), -ones(length(per_lead),1); ...
                         ], ...
                         nDofsPerBorder*nVoigt, nDofs);
            Ct =  [CtPer; CtDirPer; CtDir];
        end

        function RHSC = computeMicroDisplMonolithicRHS(obj, iVoigt, nVoigt)
            nDofs = obj.dirichletFun.nDofs;
            basis   = diag(ones(nVoigt,1));
            vstrain = basis(iVoigt,:);
            s.mesh           = obj.mesh;
            s.dirichlet_dofs = obj.dirichlet_dofs;
            s.ndofs          = nDofs;
            s.vstrain        = vstrain;
            comp = MicroDispMonolithicBCApplier(s);
            if iVoigt == 1 || iVoigt == 2
                RHSDir    = zeros(6, 1);
                RHSDirPer = -ones(2,1);
            else
                RHSDir    = zeros(4, 1);
                RHSDir(3) = 0.5;
                RHSDir(4) = 0.5;
                RHSDirPer = -ones(2,1);
            end
            zerosRHS = zeros(nDofs,1);
            nDofsPerBorder = length(obj.periodic_leader)/4; % 4 because 2D
            per_vec = zeros(nDofsPerBorder*3, 1);
            per_vec(1:nDofsPerBorder) = -vstrain(1);
            per_vec(nDofsPerBorder+1:2*nDofsPerBorder) = -vstrain(3);
            per_vec(2*nDofsPerBorder+1:end) = -vstrain(2);
            [RHSDir, RHSDirPer] = comp.getRHSVector();
            RHSC = [zerosRHS; per_vec; RHSDirPer; RHSDir];
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            inBC = cParams.boundaryConditions;
            obj.mesh = cParams.mesh;
            obj.dirichletFun    = inBC.dirichletFun;
            obj.dirichlet_dofs  = inBC.dirichlet_dofs;
            obj.dirichlet_vals  = inBC.dirichlet_vals;
            obj.pointloadFun    = inBC.pointloadFun;
            obj.pointload_dofs  = inBC.pointload_dofs;
            obj.pointload_vals  = inBC.pointload_vals;
            obj.periodic_leader = inBC.periodic_leader;
            obj.periodic_follower = inBC.periodic_follower;
        end

    end
    
end