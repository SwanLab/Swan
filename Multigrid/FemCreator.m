classdef FemCreator < handle
    
    properties (Access = public)
        LHS
        RHS
        bc
        solver
    end
    
    properties (Access = private)
        coarseMeshes
        nDimf
        nLevel
        coarseBc
        coarseDispFun
        coarseMaterial
        type
        scale
        pdim
        tol
    end
    
    methods (Access = public)
        function obj = FemCreator(cParams)
            obj.init(cParams);
            obj.createLHS();
            obj.createRHS();
            obj.createSolver();
            obj.createSolverCoarse();
        end
    end
    
    methods (Access = private)
        function init (obj,cParams)
            obj.coarseMeshes = cParams.coarseMeshes;
            obj.nDimf        = cParams.nDimf;
            obj.nLevel       = cParams.nLevel;
            obj.type         = cParams.type;
            obj.scale        = cParams.scale;
            obj.pdim         = cParams.pdim;
            obj.tol          = cParams.tol;
        end
        
        function createLHS(obj)
            for i = 1:obj.nLevel
                
                u          = P1Function.create(obj.coarseMeshes{i}, obj.nDimf);
                obj.bc{i}  = obj.createBoundaryConditions(obj.coarseMeshes{i},u);
                mat        = obj.createMaterial(obj.coarseMeshes{i});
                m          = obj.coarseMeshes{i};

                obj.LHS{i} = obj.computeKred(m,mat,u,obj.bc{i});
            end
        end
        
        function createRHS(obj)
            for i = 1:obj.nLevel
                
                u          = P1Function.create(obj.coarseMeshes{i}, obj.nDimf);
                obj.bc{i}  = obj.createBoundaryConditions(obj.coarseMeshes{i},u);
                m          = obj.coarseMeshes{i};

                obj.RHS{i} = obj.computeFred(m,u,obj.bc{i});

            end
        end
        
        function bc = createBoundaryConditions(obj,mesh,disp)
            rawBc       = obj.createRawBoundaryConditions(mesh);
            dim         = getFunDims(obj,mesh,disp);
            rawBc.ndimf = dim.ndimf;
            rawBc.ndofs = dim.ndofs;
            s.mesh      = mesh;
            s.scale     = 'MACRO';
            s.bc        = {rawBc};
            s.ndofs     = dim.ndofs;
            bc          = BoundaryConditions(s);
            bc.compute();
        end
        
        function dim = getFunDims(obj,mesh,disp)
            d.ndimf     = disp.ndimf;
            d.nnodes    = size(disp.fValues, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = mesh.nnodeElem; 
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim         = d;
        end
        
        function bc = createRawBoundaryConditions(obj,mesh)
            dirichletNodes            = abs(mesh.coord(:,1)-0) < 1e-12;
            rightSide                 = max(mesh.coord(:,1));
            isInRight                 = abs(mesh.coord(:,1)-rightSide)< 1e-12;
            forceNodes                = isInRight;
            nodes                     = 1:mesh.nnodes;
            bcDir                     = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            nodesdir                  = size(nodes(dirichletNodes),2);
            bcDir(1:nodesdir,end+1)   = 1;
            bcDir(nodesdir+1:end,end) = 2;
            bcDir(:,end+1)            = 0;
            bc.dirichlet              = bcDir;
            bc.pointload(:,1)         = nodes(forceNodes);
            bc.pointload(:,2)         = 2;
            bc.pointload(:,3)         = -1/length(forceNodes);
        end
        
        function mat = createMaterial(obj,mesh)
            s.mesh  = mesh;
            s.type  = obj.type;
            s.scale = obj.scale;
            ngaus   = 1;
            Id      = ones(mesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = obj.pdim;
            s.nelem = mesh.nelem;
            s.mesh  = mesh;
            s.kappa = .9107*Id;
            s.mu    = .3446*Id;
            mat     = Material.create(s);
            mat.compute(s);
        end
        
        function Kred = computeKred(obj,m,mat,u,bc)
            K    = obj.computeStiffnessMatrix(m,mat,u);
            Kred = bc.fullToReducedMatrix(K);
        end
        
        function k = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.fun      = displacementFun;
            s.material = material;
            lhs        = LHSintegrator.create(s);
            k          = lhs.compute();
        end
        
        function Fred = computeFred(obj,m,u,bc)
            b    = obj.computeRHS(m,u,bc);
            Fred = bc.fullToReducedVector(b);
        end
        
        function RHS = computeRHS(obj,mesh,dispFun,boundaryConditions)
            dim.ndimf     = dispFun.ndimf;
            dim.nnodes    = size(dispFun.fValues, 1);
            dim.ndofs     = dim.nnodes*dim.ndimf;
            dim.nnodeElem = mesh.nnodeElem; 
            dim.ndofsElem = dim.nnodeElem*dim.ndimf;
            c.dim         = dim;
            c.mesh        = mesh;
            c.BC          = boundaryConditions;
            RHS           = RHSintegrator_ElasticMacro(c);
            RHS           = RHS.compute();
        end
        
        function createSolver(obj)
            for i = 2:obj.nLevel+1
                s.maxIter             = 5;
                s.tol                 = obj.tol;
                s.solverType          = 'ITERATIVE';
                s.iterativeSolverType = 'CG';
                obj.solver{i}         = Solver.create(s);
            end
        end
        
        function createSolverCoarse(obj)
            s.maxIter             = 100000;
            s.tol                 = obj.tol;
            s.solverType          = 'ITERATIVE';
            s.iterativeSolverType = 'CG';
            obj.solver{1}         = Solver.create(s);
        end
        
    end
end

