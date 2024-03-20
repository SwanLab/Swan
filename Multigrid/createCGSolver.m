classdef createCGSolver < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        Property1
    end
    
    properties (Access = private)
        solverType
        iterativeSolverType
        tol
        nLevel
        LHS
        RHS
        coarseLHS
        coarseRHS
    end
    
    methods (Access = public)
        function obj = createCGSolver(cParams)
            obj.init(cParams);
            obj.createVCycle();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.solverType          = cParams.solverType;
            obj.iterativeSolverType = cParams.iterativeSolverType;
            obj.tol                 = cParams.tol;
            obj.nLevel              = cParams.nLevel;
            obj.LHS                 = cParams.LHS;
            obj.RHS                 = cParams.RHS;
            obj.coarseLHS           = cParams.coarseLHS;
            obj.coarseRHS           = cParams.coarseRHS;
            obj.coarseLHS{1,6}      = obj.LHS;
            obj.coarseRHS{1,6}      = obj.RHS;
        end
        
        function createVCycle(obj)
            u = 0*obj.coarseRHS{obj.nLevel + 1};
            level = obj.nLevel + 1;
            b = obj.RHS;
            %res_record = norm(obj.RHS(obj.nLevel),inf);
            while norm(obj.coarseRHS{obj.nLevel + 1} - obj.coarseLHS{obj.nLevel + 1} * u, inf) >= obj.tol
                u = vCycle(u,b,level);
                %res_record(i) = norm(data(nlevels).b - data(nlevels).A * u,inf);
                %i = i + 1;
            end
        end
        
        function u = vCycle(obj,u,b,level)
            if level == 1
                maxIter = 1000;
                solver = Solver.create('DIRECT');
                u = solver.solve(LHS,RHS);
              %  u = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter, level, mesh, numero,res);
            else
                
                s.maxIter = 20;
                s.solverType = obj.solverType;
                s.iterativeSolverType = obj.iterativeSolverType;
                s.LHS = obj.coarseLHS(level);
                %s.RHS = obj.coarseRHS(level);
                s.RHS = b;
                s.u = u;
                solver = Solver.create(s);
                u = solver.solve();
                %u = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter, level, mesh, numero,res);
                %r = b - data(level).A*u;
                r = b - obj.coarseLHS(leve)*u;
                [Rr] = interpolate(r,bc,data,level);
                [ur] = interpolate(u,bc,data,level);
                [er,numero,res] = vCycle(0*ur, Rr, level - 1);
                e = restriction(er,bc,data,level);
                u = u + e;
                solver = Solver.create('CG');
                u = solver.solve(LHS,RHS);                
                %[u, res, numero] = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter, level, mesh, numero,res);
             end

        end

        function fF = restriction(fC,bc,data,level)
                fF = bc{level - 1}.reducedToFullVector(fC);
                fF = reshape(fF,2,[])';
                fF = data(level - 1).T * fF; 
                fF = reshape(fF',[],1);
                fF = bc{level}.fullToReducedVector(fF);
        end

        function [fCoarse] = interpolate(fFine,bc,data,level)

                fFine = bc{level}.reducedToFullVector(fFine);
                fFine = reshape(fFine,2,[])';
                fCoarse = data(level - 1).R * fFine;
                fCoarse = reshape(fCoarse',[],1);
                fCoarse = bc{level - 1}.fullToReducedVector(fCoarse);
        end
    end
end

