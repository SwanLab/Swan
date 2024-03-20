classdef FemCreator < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        LHS
        RHS
        bc
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
    end
    
    methods (Access = public)
        function obj = FemCreator(cParams)
            obj.init(cParams)
            obj.createLHSandRHS()
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
        end
        
        function createLHSandRHS(obj)
            for i = 1:obj.nLevel
                
                u   = P1Function.create(obj.coarseMeshes{i}, obj.nDimf);
                obj.bc{i}  = obj.createBoundaryConditions(obj.coarseMeshes{i},u);
                mat = obj.createMaterial(obj.coarseMeshes{i});
                m   = obj.coarseMeshes{i};
%                 s.solverTyp             = 'ITERATIVE';
%                 s.iterativeSolverTyp    = 'CG';
%                 s.tol                   = 1e-6;
%                 s.maxIter               = 20;

                obj.LHS{i} = obj.computeStiffnessMatrix(m,mat,u);
                obj.RHS{i} = obj.createRHS(m,u,obj.bc{i});

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
        
        function k = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.fun      = displacementFun;
            s.material = material;
            lhs        = LHSintegrator.create(s);
            k          = lhs.compute();
        end
        
        function RHS = createRHS(obj,mesh,dispFun,boundaryConditions)
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
        
    end
end

