classdef BCApplier < handle
    
    % Goal: group dirichlet and neumann conditions
    % to allow multiple boundary conditions at the same time
    % Use: BCApplier.computeLinearConditionsMatrix()

    properties (Access = public)
    end
    
    properties (Access = private)
        mesh

        dirichlet_dofs, dirichlet_vals, dirichlet_domain
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
        
        function Ct = computeLinearConditionsMatrix(obj, order)
            switch order
                case 'Dirac'
                    % generalize lagrange multiplier -> dirac and stuff
                    dir_dofs = obj.dirichlet_dofs;
                    nDofs = obj.dirichletFun.nDofs;
                    nDirich = length(dir_dofs);
                    Ct = sparse(1:nDirich, dir_dofs, 1, nDirich, nDofs);
                otherwise
                    dir_dom = obj.dirichlet_domain;
                    [mesh_left2, l2g_mesh] = obj.mesh.getBoundarySubmesh(dir_dom);
        
                    dLambda = LagrangianFunction.create(mesh_left2, obj.mesh.ndim, order); % !!
                    uFun    = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1'); % !!
        
                    b.mesh  = mesh_left2;
                    b.test  = dLambda;
                    b.trial = uFun.restrictTo(dir_dom);
                    b.type  = 'MassMatrix';
                    lhs = LHSintegrator.create(b);
                    LHS = lhs.compute();
                    [iLoc,jLoc,vals] = find(LHS);

                    l2g_dof = ((l2g_mesh*dLambda.ndimf)' - ((dLambda.ndimf-1):-1:0))';
                    l2g_dof = l2g_dof(:);
                    jGlob = l2g_dof(jLoc);
                    Ct = sparse(iLoc,jGlob,vals, dLambda.nDofs, uFun.nDofs);
            end
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

        function rVec = fullToReducedVectorDirichlet(obj,fVec)
            dofs      = 1:1:obj.dirichletFun.nDofs;
            free_dofs = setdiff(dofs, obj.dirichlet_dofs);
            rVec      = fVec(free_dofs);
        end

        function rMat = fullToReducedMatrixDirichlet(obj,fMat)
            dofs      = 1:1:obj.dirichletFun.nDofs;
            free_dofs = setdiff(dofs, obj.dirichlet_dofs);
            rMat      = fMat(free_dofs,free_dofs);
        end

        function fVec = reducedToFullVectorDirichlet(obj,rVec)
            dofs                     = 1:1:obj.dirichletFun.nDofs;
            free_dofs                = setdiff(dofs, obj.dirichlet_dofs);
            fVec                     = zeros(obj.dirichletFun.nDofs,1);
            fVec(free_dofs)          = rVec;
            fVec(obj.dirichlet_dofs) = obj.dirichlet_vals;
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            inBC = cParams.boundaryConditions;
            obj.mesh = cParams.mesh;
            obj.dirichletFun    = inBC.dirichletFun;
            obj.dirichlet_dofs  = inBC.dirichlet_dofs;
            obj.dirichlet_vals  = inBC.dirichlet_vals;
            obj.dirichlet_domain = inBC.dirichlet_domain;
            obj.pointloadFun    = inBC.pointloadFun;
            obj.pointload_dofs  = inBC.pointload_dofs;
            obj.pointload_vals  = inBC.pointload_vals;
            obj.periodic_leader = inBC.periodic_leader;
            obj.periodic_follower = inBC.periodic_follower;
        end

    end
    
end