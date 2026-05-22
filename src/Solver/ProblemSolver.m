classdef ProblemSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        type, mode
        boundaryConditions
        BCApplier
        solver
    end
    
    properties (Access = private)
        C
    end
    
    methods (Access = public)
        
        function obj = ProblemSolver(cParams)
            obj.init(cParams);
            obj.computeRollerLinkage(cParams);
        end

        function [u,L] = solve(obj,cParams)
            [LHS, RHS] = obj.computeMatrices(cParams);
            sol        = obj.solver.solve(LHS, RHS);
            [u, L]     = obj.cleanupSolution(sol,cParams.stiffness);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.type               = cParams.solverType;
            obj.mode               = cParams.solverMode;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.BCApplier          = cParams.BCApplier;
            obj.solver             = cParams.solver;
        end

        function [LHS, RHS] = computeMatrices(obj,cParams)
            LHS = obj.assembleLHS(cParams);
            RHS = obj.assembleRHS(cParams);
        end

        function [u, L] = cleanupSolution(obj,sol,stiffness)
            %bcapp = obj.BCApplier;
            bcs   = obj.boundaryConditions;
            %hasPeriodic = ~isequal(bcs.periodic_leader, []);
            switch true
                case strcmp(obj.type, 'MONOLITHIC')
                    nDisp = size(stiffness,1);
                    u = sol(1:nDisp, :);
                    L = -sol( (nDisp+1):end, : );
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'DISP')
                    dofs = 1:size(stiffness,1);
                    free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
                    u = zeros(size(stiffness,1), 1);
                    u(free_dofs) = sol;
                    u(bcs.dirichlet_dofs) = bcs.dirichlet_vals;
                    L = [];
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'ROLLER')
                    dofs = 1:size(stiffness,1);
                    free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
                    u = zeros(size(stiffness,1), 1);
                    u(free_dofs) = sol(1:1:length(free_dofs));
                    u(bcs.dirichlet_dofs) = bcs.dirichlet_vals;
                    L = [];
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'FLUC')
                    lead = bcs.periodic_leader;
                    fllw = bcs.periodic_follower;
                    drch = bcs.dirichlet_dofs;
                    dofs = 1:size(stiffness,1);
                    free = setdiff(dofs, [lead; fllw; drch]);
                    u = zeros(length(dofs),1);
                    u(free) = sol(1:1:size(free,2));
                    u(lead) = sol(size(free,2)+1:1:size(sol,1));
                    u(fllw) = u(lead);
                    u(drch) = bcs.dirichlet_vals;
                    L = [];
                otherwise
                    u = [];
                    L = [];
            end

            % maybe return a p1function or whatever
        end

        function LHS = assembleLHS(obj,cParams)
            stiffness = cParams.stiffness;
            bcapp = obj.BCApplier;
            bcs   = obj.boundaryConditions;
            hasPeriodic = ~isequal(bcs.periodic_leader, []);

            switch true
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'DISP')
                    if ~hasPeriodic
                        % Ct = bcapp.computeLinearConditionsMatrix('P1');
                        Ct = bcapp.computeLinearConditionsMatrix('Dirac');
                        C   = Ct';
                        nC  = size(Ct,1);
                        Z   = zeros(nC);
                        Km  = stiffness;
                        LHS = [Km C; C' Z];
                    else
                        % Micro
                        iV = cParams.iBase;
                        nV = cParams.nBasis;
                        % CtDir = bcapp.computeLinearConditionsMatrix();
                        % CtPer = bcapp.computeLinearPeriodicConditionsMatrix();
                        Ct = bcapp.computeSingleDirichletPeriodicCondition(iV, nV);
                        C   = Ct';
                        nC  = size(Ct,1);
                        Z   = zeros(nC);
                        Km  = stiffness;
                        LHS = [Km C; C' Z];
                    end
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'DISP')
                    dofs = 1:size(stiffness,1);
                    free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
                    LHS = stiffness(free_dofs, free_dofs);
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'ROLLER')
                    dofs = 1:size(stiffness,1);
                    free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
                    Kff = stiffness(free_dofs, free_dofs);
                    Cf = obj.C(free_dofs,:);
                    Z = sparse(size(Cf,2),size(Cf,2));
                    LHS = [Kff, Cf; Cf', Z];
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'FLUC')
                    CtDir = bcapp.computeLinearConditionsMatrix('Dirac');
                    CtPer = bcapp.computeLinearPeriodicConditionsMatrix();
                    Ct = [CtPer; CtDir];
                    C   = Ct';
                    nC  = size(Ct,1);
                    Z   = zeros(nC);
                    Km  = stiffness;
                    LHS = [Km C; C' Z];
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'FLUC')
                    lead = bcs.periodic_leader;
                    fllw = bcs.periodic_follower;
                    drch = bcs.dirichlet_dofs;
                    dofs = 1:size(stiffness,1);
                    free = setdiff(dofs, [lead; fllw; drch]);
                    A = stiffness;
                    A_II = A(free,free);
                    A_IP = A(free,lead) + A(free,fllw); %Grouping P and Q nodal values
                    A_PI = A(lead,free) + A(fllw,free); % Adding P  and Q equation
                    A_PP = A(lead,lead) + A(lead,fllw) + A(fllw,lead) + A(fllw,fllw); % Adding and grouping
                    LHS = [A_II, A_IP; A_PI, A_PP];
            end

        end

        function RHS = assembleRHS(obj,cParams)
            stiffness = cParams.stiffness;
            forces    = cParams.forces;
            bcapp = obj.BCApplier;
            bcs   = obj.boundaryConditions;
            hasPeriodic = ~isequal(bcs.periodic_leader, []);

            switch true
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'DISP')
                    if ~hasPeriodic
                        nCases = size(forces,2);
                        % lambda = zeros(size(obj.lhs,1) - size(forces,1), 1);
                        lambda = bcs.dirichlet_vals;
                        Ct = repmat(lambda, [1 nCases]);
                        RHS = [forces; Ct];
                    else
                        iV = cParams.iBase;
                        nV = cParams.nBasis;
                        RHS = bcapp.computeMicroDisplMonolithicRHS(iV, nV);
                    end
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'DISP')
                    dofs = 1:size(stiffness,1);
                    free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
                    RHS = forces(free_dofs);
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'ROLLER')
                    dofs = 1:size(stiffness,1);
                    free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
                    Ff = forces(free_dofs);
                    Z = sparse(size(obj.C,2),1);
                    RHS = [Ff;Z];
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'FLUC')
                    nPer = length(bcs.periodic_leader);
                    RHS = [forces; zeros(nPer,1); bcs.dirichlet_vals];
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'FLUC')
                    lead = bcs.periodic_leader;
                    fllw = bcs.periodic_follower;
                    drch = bcs.dirichlet_dofs;
                    dofs = 1:size(stiffness,1);
                    free = setdiff(dofs, [lead; fllw; drch]);
                    b = forces;
                    b_I = b(free);
                    b_P = b(lead) + b(fllw);
                    RHS = [b_I; b_P];
            end

        end

        function C = computeRollerLinkage(obj,s)
            if isfield(s,'rollerMesh')
                bDirMesh = s.rollerMesh;

                n0 = squeeze(bDirMesh.getNormals())';

                nP0 = LagrangianFunction.create(bDirMesh,3,'P0');
                nP0.setFValues(n0);

                sF.trial = LagrangianFunction.create(bDirMesh,3,'P1');
                sF.mesh = bDirMesh;
                filter = FilterLump(sF);
                n0Reg = filter.compute(nP0,2);
                n0 = n0Reg.fValues;
                n0 = n0./sqrt(n0(:,1).^2+n0(:,2).^2+n0(:,3).^2);
                n0 = n0';

                jC = [1:n0Reg.nDofs/n0Reg.ndimf;1:n0Reg.nDofs/n0Reg.ndimf;1:n0Reg.nDofs/n0Reg.ndimf];
                jC = reshape(jC,[],1);
                iC = 1:n0Reg.nDofs;
                iC = iC';
                CLoc = sparse(iC,jC,n0(:),n0Reg.nDofs,n0Reg.nDofs/n0Reg.ndimf);

                obj.C = sparse(s.mesh.nnodes*3,n0Reg.nDofs/n0Reg.ndimf);
                iC = s.rollerNodes;
                iC = repmat(iC,[1,3])';
                iC(1,:) = 3*iC(1,:)-2;
                iC(2,:) = 3*iC(2,:)-1;
                iC(3,:) = 3*iC(3,:);
                iC = iC(:);
                obj.C(iC,unique(jC)) = CLoc;
            end
        end

    end

end