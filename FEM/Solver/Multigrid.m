classdef Multigrid < handle


    properties (Access = private)
        nDimf
        data
        nLevel
        mesh
        material
        boundaryConditions
        tol
        type
        scale  
        pdim
        RHS
        LHS
        coarseMeshes
        interpolator
        coarseBc
        coarseLHS
        coarseRHS
        solverType
        iterativeSolverType
    end

    methods (Access = public)

        function obj = Multigrid(cParams)
            obj.init(cParams);            
            obj.createFEMlevel();
        end

        function u = solve(obj)
            obj.coarseBc{1,obj.nLevel + 1}  = obj.boundaryConditions;
            obj.coarseLHS{1,obj.nLevel + 1} = obj.LHS;
            obj.coarseRHS{1,obj.nLevel + 1} = obj.RHS;
            u                               = 0*obj.coarseRHS{obj.nLevel + 1};
            level                           = obj.nLevel + 1;
            b                               = obj.RHS;
            
            while norm(obj.coarseRHS{obj.nLevel + 1} - obj.coarseLHS{obj.nLevel + 1} * u, inf) >= obj.tol
                u = obj.vCycle(u,b,level);
            end

        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.coarseMeshes       = cParams.coarseMeshes;
            obj.interpolator       = cParams.interpolator;
            obj.tol                = cParams.tol;
            obj.nLevel             = cParams.nLevel;
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.bc;
            obj.material           = cParams.material;
            obj.LHS                = cParams.LHS;
            obj.RHS                = cParams.RHS;
            obj.type               = cParams.type ;%'ELASTIC';
            obj.scale              = cParams.scale; %'MACRO';
            obj.pdim               = cParams.dim; %'2D';
            obj.nDimf              = cParams.nDimf;
        end

        function createFEMlevel(obj)
            s.coarseMeshes = obj.coarseMeshes;
            s.nDimf        = obj.nDimf;
            s.nLevel       = obj.nLevel;
            s.type         = obj.type;
            s.scale        = obj.scale;
            s.pdim         = obj.pdim;
            FEM            = FemCreator(s);
            obj.coarseBc   = FEM.bc;
            obj.coarseLHS  = FEM.LHS;
            obj.coarseRHS  = FEM.RHS;
        end
        
        function u = vCycle(obj,u,b,level)
            if level == 1
                maxIter               = 100000;
                LHS                   = obj.coarseLHS{level};
                RHS                   = b;
                s.solverType          = 'ITERATIVE';
                s.iterativeSolverType = 'CG';
                solver                = Solver.create(s);
                u                     = solver.solve(LHS,RHS,maxIter,u);
            else
                maxIter = 20;
                LHS                   = obj.coarseLHS{level};
                RHS                   = b;
                s.solverType          = 'ITERATIVE';
                s.iterativeSolverType = 'CG';
                solver                = Solver.create(s);
                u                     = solver.solve(LHS,RHS,maxIter,u);
                r                     = b - obj.coarseLHS{level}*u;
                Rr                    = obj.interpolate(r,obj.coarseBc,obj.interpolator,level);
                ur]+                  = obj.interpolate(u,obj.coarseBc,obj.interpolator,level);
                er                    = obj.vCycle(0*ur, Rr, level - 1);
                e                     = obj.restriction(er,obj.coarseBc,obj.interpolator,level);
                u                     = u + e;
                maxIter               = 20;
                LHS                   = obj.coarseLHS{level};
                RHS                   = b;
                s.solverType          = 'ITERATIVE';
                s.iterativeSolverType = 'CG';
                solver                = Solver.create(s);
                u                     = solver.solve(LHS,RHS,maxIter,u);                
             end

        end

        function fF = restriction(obj,fC,bc,I,level)
                fF = bc{level - 1}.reducedToFullVector(fC);
                fF = reshape(fF,2,[])';
                fF = I{level - 1} * fF; 
                fF = reshape(fF',[],1);
                fF = bc{level}.fullToReducedVector(fF);
        end

        function fCoarse = interpolate(obj,fFine,bc,I,level)
                fFine   = bc{level}.reducedToFullVector(fFine);
                fFine   = reshape(fFine,2,[])';
                fCoarse = I{level - 1}' * fFine;
                fCoarse = reshape(fCoarse',[],1);
                fCoarse = bc{level - 1}.fullToReducedVector(fCoarse);
        end
        
    end
end

