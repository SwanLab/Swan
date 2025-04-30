classdef ProblemSolverOld < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        type, mode
        boundaryConditions
        BCApplier
        solver
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ProblemSolver(cParams)
            obj.init(cParams);
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
                        iV = cParams.iVoigt;
                        nV = cParams.nVoigt;
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
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'FLUC')
                    CtDir = bcapp.computeLinearConditionsMatrix('Dirac');
                    CtPer = bcapp.computeLinearPeriodicConditionsMatrix();
                    CtCs  = bcapp.computeLinearConditionsMatrix('Analytical');
                    Cv    = bcapp.computeVoluMatrix();
                    Ct      = [CtPer; CtDir; CtCs';Cv'];
                    C       = Ct';
                    nC      = size(Ct,1);
                    Z       = zeros(nC);
                    Km      = stiffness;
                    LHS     = [Km C; C' Z];

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
                        iV = cParams.iVoigt;
                        nV = cParams.nVoigt;
                        RHS = bcapp.computeMicroDisplMonolithicRHS(iV, nV);
                    end
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'DISP')
                    dofs = 1:size(stiffness,1);
                    free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
                    RHS = forces(free_dofs);
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'FLUC')
                    nPer = length(bcs.periodic_leader);
                    RHS = [forces; zeros(nPer,1); bcs.dirichlet_vals; [0;0;0;0;0]];
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
        
    end
    
end