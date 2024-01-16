classdef ProblemSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        type, mode
        stiffness
        forces
        boundaryConditions
        BCApplier
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ProblemSolver(cParams)
            obj.init(cParams)
        end

        function u = solve(obj)
            [LHS, RHS] = obj.computeMatrices();
            sol        = obj.solveSystem(LHS, RHS);
            [u, L]     = obj.cleanupSolution(sol);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.type               = cParams.solverType;
            obj.mode               = cParams.solverMode;
            obj.stiffness          = cParams.stiffness;
            obj.forces             = cParams.forces;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.BCApplier          = cParams.BCApplier;
        end

        function [LHS, RHS] = computeMatrices(obj)
            LHS  = obj.assembleLHS();
            RHS  = obj.assembleRHS();
        end

        function sol = solveSystem(obj, LHS, RHS)
            sol = LHS\RHS;
        end

        function [u, L] = cleanupSolution(obj,sol)
            bcapp = obj.BCApplier;
            hasPeriodic = ~isequal(bcapp.periodic_leader, []);
            switch true
                case strcmp(obj.type, 'MONOLITHIC')
                    nDisp = size(obj.stiffness,1);
                    u = sol(1:nDisp, :);
                    L = -sol( (nDisp+1):end, : );
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'DISP')
                    dofs = 1:size(obj.stiffness);
                    free_dofs = setdiff(dofs, bcapp.dirichlet_dofs);
                    u = zeros(size(obj.stiffness,1), 1);
                    u(free_dofs) = sol;
                    u(bcapp.dirichlet_dofs) = bcapp.dirichlet_vals;
                    L = [];
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'FLUC')
                    lead = bcapp.periodic_leader;
                    fllw = bcapp.periodic_follower;
                    drch = bcapp.dirichlet_dofs;
                    dofs = 1:size(obj.stiffness);
                    free = setdiff(dofs, [lead; fllw; drch]);
                    u = zeros(length(dofs),1);
                    u(free) = sol(1:1:size(free,2));
                    u(lead) = sol(size(free,2)+1:1:size(sol,1));
                    u(fllw) = u(lead);
                    L = [];
                otherwise
                    u = [];
                    L = [];
            end

            % maybe return a p1function or whatever
        end

        function LHS = assembleLHS(obj)
            bcapp = obj.BCApplier;
            hasPeriodic = ~isequal(bcapp.periodic_leader, []);

            switch true
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'DISP')
                    if ~hasPeriodic
                        Ct = bcapp.computeLinearConditionsMatrix();
                        C   = Ct';
                        nC  = size(Ct,1);
                        Z   = zeros(nC);
                        Km  = obj.stiffness;
                        LHS = [Km C; C' Z];
                    else
                        CtDir = bcapp.computeLinearConditionsMatrix();
                        CtPer = bcapp.computeLinearPeriodicConditionsMatrix();
                        CtPerDir = bcapp.computeSingleDirichletPeriodicCondition();
                    end
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'DISP')
                    dofs = 1:size(obj.stiffness);
                    free_dofs = setdiff(dofs, bcapp.dirichlet_dofs);
                    LHS = obj.stiffness(free_dofs, free_dofs);
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'FLUC')
                    CtDir = bcapp.computeLinearConditionsMatrix();
                    CtPer = bcapp.computeLinearPeriodicConditionsMatrix();
                    Ct = [CtPer; CtDir];
                    C   = Ct';
                    nC  = size(Ct,1);
                    Z   = zeros(nC);
                    Km  = obj.stiffness;
                    LHS = [Km C; C' Z];
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'FLUC')
                    lead = bcapp.periodic_leader;
                    fllw = bcapp.periodic_follower;
                    drch = bcapp.dirichlet_dofs;
                    dofs = 1:size(obj.stiffness);
                    free = setdiff(dofs, [lead; fllw; drch]);
                    A = obj.stiffness;
                    A_II = A(free,free);
                    A_IP = A(free,lead) + A(free,fllw); %Grouping P and Q nodal values
                    A_PI = A(lead,free) + A(fllw,free); % Adding P  and Q equation
                    A_PP = A(lead,lead) + A(lead,fllw) + A(fllw,lead) + A(fllw,fllw); % Adding and grouping      
                    LHS = [A_II, A_IP; A_PI, A_PP];
            end

        end

        function RHS = assembleRHS(obj)
            bcapp = obj.BCApplier;
            hasPeriodic = ~isequal(bcapp.periodic_leader, []);

            switch true
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'DISP')
                    if ~hasPeriodic
                        dir_vals = bcapp.dirichlet_vals;
                        nCases = size(obj.forces,2);
                        Ct = repmat(dir_vals, [1 nCases]);
                        RHS = [obj.forces; Ct];
                    else

                    end
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'DISP')
                    dofs = 1:size(obj.stiffness);
                    free_dofs = setdiff(dofs, bcapp.dirichlet_dofs);
                    RHS = obj.forces(free_dofs);
                case strcmp(obj.type, 'MONOLITHIC') && strcmp(obj.mode, 'FLUC')
                    nPer = length(bcapp.periodic_leader);
                    RHS = [obj.forces; zeros(nPer,1); bcapp.dirichlet_vals];
                case strcmp(obj.type, 'REDUCED') && strcmp(obj.mode, 'FLUC')
                    lead = bcapp.periodic_leader;
                    fllw = bcapp.periodic_follower;
                    drch = bcapp.dirichlet_dofs;
                    dofs = 1:size(obj.stiffness);
                    free = setdiff(dofs, [lead; fllw; drch]);
                    b = obj.forces;
                    b_I = b(free);
                    b_P = b(lead) + b(fllw);
                    RHS = [b_I; b_P];
            end

        end
        
    end
    
end