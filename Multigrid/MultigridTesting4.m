classdef MultigridTesting4 < handle
    
    properties (Access = public)
        u
    end
    
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
            tic
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
            s.RHS                 = obj.computeFred(mF,s.dispFun,s.bc);
            s.LHS                 = obj.computeKred(mF,s.material,s.dispFun,s.bc);
            s.type                = 'ELASTIC';
            s.scale               = 'MACRO';
            s.dim                 = '3D';
            s.solverType          = 'ITERATIVE';
            s.iterativeSolverType = 'MULTIGRID';
            s.tol                 = 1e-6;
            s.nLevel              = obj.nLevel;
            s.nDimf               = obj.nDimf;
            solver                = Solver.create(s);
            obj.u                 = solver.solve();
            toc
            
            obj.postProcess();
            obj.plotRes(obj.u,mF,s.bc);
        end
    end

    methods (Access = private)
        
        function init(obj)
            obj.nDimf        = 3;
            obj.nbasis       = 20;
            obj.functionType = 'P1';
            obj.nLevel       = 3;
        end

        function createMultiLevelMesh(obj)
           s.nX               = 1;
           s.nY               = 1;
           s.nZ               = 1;
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
            s.pdim  = '3D';
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
        
        function LHS = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.fun      = displacementFun;
            s.material = material;
            lhs        = LHSintegrator.create(s);
            LHS        = lhs.compute();
        end 
        
        function Fred = computeFred(obj,m,u,bc)
            b  = obj.createRHS(m,u,bc);
            Fred = bc.fullToReducedVector(b);
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
        
        function postProcess(obj)
            DOFs = [1056 9312 25760 50400 83232];
            time2 = [0.051821 0.072445 0.113856 0.229354 0.282576];
            time5 = [0.00513 0.020559 0.046073 0.090326 0.130545];
            time10 = [0.004523 0.031546 0.067093 0.135769 0.16341];
            time20 = [0.005 0.046528 0.084207 0.174044 0.292117];
            
            figure(1)
            plot(DOFs,time2)
            hold on
            plot(DOFs,time5)
            plot(DOFs,time10)
            plot(DOFs,time20)
            title('Time vs DOFs')
            xlabel('DOFs')
            ylabel('Time(s)')
            
            legend('2 iters per level','5 iters per level','10 iters per level','20 iters per level','location','northwest')

            time2Total = [0.887724, 1.18978, 7.160798, 38.034618, 155.408853];
            time5Total = [0.33394, 1.047098, 7.062304, 39.194001, 163.584521];
            time10Total = [0.259715, 1.041644, 7.06799, 39.915648, 159.657016];
            time20Total = [0.266481, 0.97331, 7.010966, 38.85078, 151.357202];

            figure(2)
            plot(DOFs,time2Total)
            hold on
            plot(DOFs,time5Total)
            plot(DOFs,time10Total)
            plot(DOFs,time20Total)
            title('Total Time vs DOFs')
            xlabel('DOFs')
            ylabel('Time(s)')
            
            legend('2 iters per level','5 iters per level','10 iters per level','20 iters per level','location','northwest')
        end

        function plotRes(obj,res,mesh,bc,numItr)
            xFull = bc.reducedToFullVector(res);
            s.fValues = reshape(xFull,3,[])';
            s.mesh = mesh;
            %s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print('uPrueva','Paraview')
            fclose('all');
        end
    end

end