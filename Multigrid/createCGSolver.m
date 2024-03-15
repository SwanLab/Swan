classdef createCGSolver
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
            obj.init(cParams)
            obj.createVCycle();
        end
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.solverType = cParams.solverType;
            obj.iterativeSolverType = cParams.iterativeSolverType;
            obj.tol = cParams.tol;
            obj.nLevel = cParams.nLevel;
            obj.LHS = cParams.LHS;
            obj.RHS = cParams.RHS;
            obj.coarseLHS = cParams.coarseLHS;
            obj.coarseRHS = cParams.coarseRHS;
        end
        
        function createVCycle(obj)
            u = 0*obj.RHS(obj.nLevel);
            res_record = norm(obj.RHS(obj.nLevel),inf);
            while norm(obj.RHS(obj.nLevel) - obj.LHS(obj.nLevel) * u, inf) >= obj.tol
                [u,numero] = vcycle(u, data(nlevels).b, data, vdown, vup, nlevels, bc, numero, mesh, meshType);
                res_record(i) = norm(data(nlevels).b - data(nlevels).A * u,inf);
                i = i + 1;
            end
        end
        
        function vCycle()
            if level == 1
                maxIter = 1000;
                [u, res, numero] = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter, level, mesh, numero,res);
            else 
                maxIter = 20;
                [u, res, numero] = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter, level, mesh, numero,res);
                r = b - data(level).A*u;
                [Rr] = interpolate(r,bc,data,level);
                [ur] = interpolate(u,bc,data,level);
                [er,numero,res] = vcycle(0*ur, Rr, data, vdown, vup, level - 1, bc, numero, mesh, meshType,res);
                e = restriction(er,bc,data,level);
                u = u + e;
                [u, res, numero] = conjugateGradient_Solver(data(level).A,b,u,meshType,maxIter, level, mesh, numero,res);
             end

        end

        function fFine = restriction(fCoarse,bc,data,level)
                fFine = bc{level - 1}.reducedToFullVector(fCoarse);
                fFine = reshape(fFine,2,[])';
                fFine = data(level - 1).T * fFine; 
                fFine = reshape(fFine',[],1);
                fFine = bc{level}.fullToReducedVector(fFine);
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

