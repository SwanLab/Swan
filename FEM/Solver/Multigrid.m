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
        coarseLHS
        coarseRHS
        
%         solver
    end

    methods (Access = public)

        function obj = Multigrid(cParams)
            obj.init(cParams);
%             s.type                = 'ELASTIC';
%             s.scale               = 'MACRO';
%             s.dim                 = '2D';
%             s.solverTyp           = 'ITERATIVE';
%             s.iterativeSolverType = 'CG';
            
            obj.createFEMlevel();
            obj.createSolver();


            obj.createData();

        end

        function u = solve(obj)
            while norm(data(nlevels).b - data(nlevels).A * u, inf) >= obj.tol
                [u,numero] = vcycle(u, data(nlevels).b, data, vdown, vup, nlevels, bc, numero, mesh, meshType);
                res_record(i) = norm(data(nlevels).b - data(nlevels).A * u,inf);
                i = i + 1;
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
%             s.solverType           = 'DIRECT';
%             obj.solver             = Solver.create(s);
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
            obj.coarseLHS  = FEM.LHS;
            obj.coarseRHS  = FEM.RHS;
        end
        
        function createSolver(obj)
            s.solverType          = 'ITERATIVE';
            s.iterativeSolverType = 'CG';
            s.tol                 = obj.tol;
            s.nLevel              = obj.nLevel;
            s.LHS                 = obj.LHS;
            s.RHS                 = obj.RHS;
            s.coarseLHS           = obj.coarseLHS;
            s.coarseRHS           = obj.coarseRHS;
            CGSolverCreator(s);
            
        end
    end
end

