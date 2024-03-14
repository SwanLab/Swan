classdef MultigridTesting4 < handle
    
    
    properties (Access = private)
        nLevel        
        multiLevelMesh


        nDimf
        boundaryConditionsFine
        fineMaterial
        quad
        
        functionType
        nbasis
        
        Lt
        FineK
        
        FineKred
        coarseMesh
        CoarseK
        CoarseKred
        coarseMaterial
        boundaryConditionsCoarse
        fineMeshCoord
        fineMeshConnec
        fineMesh
        FineDispFun
        CoarseDispFun
        FineFred
        CoarseFred
        data
        
        
        I
        material
        boundaryConditions
        dispFun
        Kred
        Fred
        fem
    end
    
    methods (Access = public)
        
        function obj = MultigridTesting4()
            close all;
            addpath(genpath(fileparts(mfilename('fullpath'))))
            obj.init();
            obj.createMultiLevelMesh();
            mF           = obj.multiLevelMesh.mesh;
            coarseMeshes = obj.multiLevelMesh.coarseMeshes;
            interpolator = obj.multiLevelMesh.interpolator;
             

            s.mesh                = mF;
            s.coarseMeshes        = coarseMeshes;
            s.interpolator        = interpolator;  
            s.bc                  = obj.createBoundaryConditions(mF);
            s.material            = obj.createMaterial(mF);
            s.dispFun             = P1Function.create(mF, obj.nDimf);
            s.RHS                 = obj.createRHS(mF,s.dispFun,s.bc);
            s.LHS                 = obj.computeStiffnessMatrix(mF,s.material,s.dispFun);
            s.type                = 'ELASTIC';
            s.scale               = 'MACRO';
            s.dim                 = '2D';
            s.solverType          = 'ITERATIVE';
            s.iterativeSolverType = 'MULTIGRID';
            s.tol                 = 1e-6;
            s.nLevel              = 5;
            s.nDimf               = 2;
            solver                = Solver.create(s);
            u                     = solver.solve();   
        end
    end

    methods (Access = private)
        
        function init(obj)
            obj.nDimf        = 2;
            obj.nbasis       = 20;
            obj.functionType = 'P1';
            obj.nLevel       = 5;
        end

        function createMultiLevelMesh(obj)
           s.nX               = 2;
           s.nY               = 2;
           s.nLevel           = obj.nLevel;
           m                  = MultilevelMesh(s);
           obj.multiLevelMesh = m;
        end
        
        function bc = createBoundaryConditions(obj,mesh)
            rawBc       = obj.createRawBoundaryConditions(mesh);
            dim         = obj.getFunDims(mesh);
            rawBc.ndimf = dim.ndimf;
            rawBc.ndofs = dim.ndofs;
            s.mesh      = obj.fineMesh;
            s.scale     = 'MACRO';
            s.bc        = {rawBc};
            s.ndofs     = dim.ndofs;
            bc          = BoundaryConditions(s);
            bc.compute();
        end
        
        function dim = getFunDims(obj,mesh)
            s.fValues   = mesh.coord;
            s.mesh      = mesh;
            disp        = P1Function(s);
            d.ndimf     = disp.ndimf;
            d.nnodes    = size(disp.fValues, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = mesh.nnodeElem; 
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim         = d;
        end
        
        function bc = createRawBoundaryConditions(obj,mesh)
            dirichletNodes = abs(mesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(mesh.coord(:,1));
            isInRight = abs(mesh.coord(:,1)-rightSide)< 1e-12;
            forceNodes = isInRight;
            nodes = 1:mesh.nnodes;
            bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            nodesdir=size(nodes(dirichletNodes),2);
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
            s.type  = 'ELASTIC';
            s.scale = 'MACRO';
            ngaus   = 1;
            Id      = ones(mesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = mesh.nelem;
            s.mesh  = mesh;
            s.kappa = .9107*Id;
            s.mu    = .3446*Id;
            mat     = Material.create(s);
            mat.compute(s);
        end

        function LHS = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.fun      = displacementFun;
            s.material = material;
            lhs        = LHSintegrator.create(s);
            LHS        = lhs.compute();
        end 
        
        function RHS = createRHS(obj,mesh,dispFun,boundaryConditions)
            dim.ndimf  = dispFun.ndimf;
            dim.nnodes = size(dispFun.fValues, 1);
            dim.ndofs  = dim.nnodes*dim.ndimf;
            dim.nnodeElem = mesh.nnodeElem; 
            dim.ndofsElem = dim.nnodeElem*dim.ndimf;
            c.dim=dim;
            c.mesh=mesh;
            c.BC = boundaryConditions;
            RHS    = RHSintegrator_ElasticMacro(c);
            RHS = RHS.compute();
        end
        
%         function createData(obj)
%             
%             for i = 1:obj.nMesh
%                 obj.data(i).p = obj.multiLevelMesh{i}.coord;
%                 obj.data(i).t = obj.multiLevelMesh{i}.connec;
%                 if i < obj.nMesh
%                     obj.data(i).T = obj.I{i};
%                     obj.data(i).R = obj.I{i}';
%                 end
%                 obj.data(i).A = obj.Kred{i};
%                 obj.data(i).b = obj.Fred{i};
%             end
%             
%         end

    end

end