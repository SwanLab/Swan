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
            mF = obj.multiLevelMesh.mesh;
             

            s.mesh = mF;
            s.bc   = obj.createBoundaryConditions(mF);
            %             obj.createBoundaryConditions(0)
            %             obj.createMaterial(0)
            s.material = obj.createMaterial(mF);
            s.dispFun  = P1Function.create(mF, obj.nDimf);
            s.RHS      = obj.createRHS(mF,s.dispFun,s.bc);
            s.LHS      = obj.computeStiffnessMatrix(mF,s.material,s.dispFun);
            s.type     = 'ELASTIC';
            s.scale    = 'MACRO';
            s.dim      = '2D';
            s.solverType           = 'ITERATIVE';
            s.iterativeSolverType = 'MULTIGRID';
            s.tol                 = 1e-6;
            s.nLevel              = 5;
            s.nDimf = 2;
            solver = Solver.create(s);
            u = solver.solve();
%             obj.fem{1} = FEM.create(s);
            %             obj.computeKred(0)
            %             obj.computeFred(0)

     

        end

        function r = getdata(obj)
            r = obj.data;
        end

        function r = getBC(obj)
            r = obj.boundaryConditions;
        end

        function r = getMesh(obj)
            r = obj.multiLevelMesh;
        end
    end

    methods (Access = private)
        
        function init(obj)
            obj.nDimf = 2;
            obj.nbasis = 20;
            obj.functionType = 'P1';
            obj.nLevel = 5;
        end

 
        

        function createMultiLevelMesh(obj)
           s.nX = 2;
           s.nY = 2;
           s.nLevel = obj.nLevel;
           m = MultilevelMesh(s);
           obj.multiLevelMesh = m;
        end

        function createMatrixInterpolation(obj,nMesh)
            mesh = obj.fem{nMesh}.getMesh();
            p = mesh.coord;
            t = mesh.connec;

            n = size(p,1);
            q = size(t,1);
            T = sparse(eye(n,n)); 
            tnew = []; j = 1;
            p_ori = p;
            for i = 1:q % this will add all the midpoints into p
                tcurr = t(i,:);
                pmid = [(p(tcurr(1),:) + p(tcurr(2),:)) / 2;
                        (p(tcurr(2),:) + p(tcurr(3),:)) / 2;
                        (p(tcurr(3),:) + p(tcurr(1),:)) / 2];
                p = [p; pmid];
            end
            
            [~,ia] = unique(p,'rows','stable');
            Ia = ia(n+1:end);
            ias = ia(n+1:end) - n ;   
            potential_tri = ceil(ias./3);
            d = 1;
            midpt_curr = [];
            
            for i = potential_tri' % now need to loop thru ia and find the triangle that 
                %corresponds to this midpoint
                tcurr = t(i,:);
                midpt_curr(1,:) = p(Ia(d),:);
                
                pmid = [(p(tcurr(1),:) + p(tcurr(2),:)) / 2;
                        (p(tcurr(2),:) + p(tcurr(3),:)) / 2;
                        (p(tcurr(3),:) + p(tcurr(1),:)) / 2];
                    
                if midpt_curr(1,:) == pmid(1,:)
                    T(n + 1, [tcurr(1),tcurr(2)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(2,:)
                    T(n + 1, [tcurr(2),tcurr(3)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(3,:)
                    T(n + 1, [tcurr(1),tcurr(3)]) = 1/2;
                end
                n = n + 1;
                d = d + 1;
            end
            obj.I{nMesh} = T;
        end
        
        function meshfine = createMesh(obj,i)
            meshcoarse = obj.fem{i}.getMesh();
            meshCoord = obj.I{i} * meshcoarse.coord;
            meshConnec = delaunayn(meshCoord);
            s.coord = meshCoord;
            s.connec = meshConnec;
%             obj.mesh{i+1} = Mesh(s);
            meshfine = Mesh(s);
        end
        
        function bc = createBoundaryConditions(obj,mesh)
            rawBc    = obj.createRawBoundaryConditions(mesh);
            dim = obj.getFunDims(mesh);
            rawBc.ndimf = dim.ndimf;
            rawBc.ndofs = dim.ndofs;
            s.mesh  = obj.fineMesh;
            s.scale = 'MACRO';
            s.bc    = {rawBc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            %obj.boundaryConditions{i+1} = bc;
        end
        
        function dim = getFunDims(obj,mesh)
            s.fValues = mesh.coord;
            s.mesh = mesh;
            disp = P1Function(s);
            d.ndimf  = disp.ndimf;
            d.nnodes = size(disp.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end
        
        function bc = createRawBoundaryConditions(obj,mesh)
            dirichletNodes = abs(mesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(mesh.coord(:,1));
            isInRight = abs(mesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(mesh.coord(:,2)-1.5) < 0.1;
            %forceNodes = isInRight & isInMiddleEdge;
            forceNodes = isInRight;
            nodes = 1:mesh.nnodes;
            bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            nodesdir=size(nodes(dirichletNodes),2);
            bcDir(1:nodesdir,end+1) = 1;
            bcDir(nodesdir+1:end,end) = 2;
            bcDir(:,end+1)=0;
            bc.dirichlet = bcDir;
            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = -1/length(forceNodes);
        end
        
        function mat = createMaterial(obj,mesh)
            s.mesh = mesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            ngaus = 1;
            Id = ones(mesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = mesh.nelem;
            s.mesh  = mesh;
            s.kappa = .9107*Id;
            s.mu    = .3446*Id;
            mat = Material.create(s);
            mat.compute(s);
%             obj.material{i+1} = mat;
        end
        
%         function LHS = computeKred(obj,i)
%             obj.dispFun{i+1} = P1Function.create(obj.mesh{i+1}, obj.nDimf);
%             LHS    = obj.computeStiffnessMatrix(obj.mesh{i+1},obj.material{i+1},obj.dispFun{i+1});
%             %obj.Kred{i+1} = obj.boundaryConditions{i+1}.fullToReducedMatrix(K);
%         end

        function LHS = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.fun      = displacementFun;
            % s.test      = displacementFun;
            % s.trial      = displacementFun;
            s.material = material;
            lhs = LHSintegrator.create(s);
            LHS   = lhs.compute();
        end 
        
%         function RHS = computeRHS(obj,i)
%             RHS  = obj.createRHS(obj.mesh{i+1},obj.dispFun{i+1},obj.boundaryConditions{i+1});
%             RHS = RHS.compute();
%             %obj.Fred{i+1} = obj.boundaryConditions{i+1}.fullToReducedVector(Fext);
%         end
        
        function RHS = createRHS(obj,mesh,dispFun,boundaryConditions)
            dim.ndimf  = dispFun.ndimf;
            dim.nnodes = size(dispFun.fValues, 1);
            dim.ndofs  = dim.nnodes*dim.ndimf;
            dim.nnodeElem = mesh.nnodeElem; % should come from interp..
            dim.ndofsElem = dim.nnodeElem*dim.ndimf;
            c.dim=dim;
            c.mesh=mesh;
            c.BC = boundaryConditions;
            RHS    = RHSintegrator_ElasticMacro(c);
            RHS = RHS.compute();
        end

        function fem = createFEM(obj)
            s.mesh     = obj.meshDomain;
            s.bc       = obj.boundaryConditions;
            s.material = obj.material;
            s.type     = 'ELASTIC';
            s.scale    = 'MACRO';
            s.dim      = '2D';
            s.solverTyp = 'ITERATIVE';
            s.iterativeSolverTyp = 'PCG';
            s.preconditionerType = 'EIFEM';
            s.tol = 1e-6;

            fem        = FEM.create(s);
%             fem.solve();
        end
        
        function createData(obj)
            
            for i = 1:obj.nMesh
                obj.data(i).p = obj.multiLevelMesh{i}.coord;
                obj.data(i).t = obj.multiLevelMesh{i}.connec;
                if i < obj.nMesh
                    obj.data(i).T = obj.I{i};
                    obj.data(i).R = obj.I{i}';
                end
                obj.data(i).A = obj.Kred{i};
                obj.data(i).b = obj.Fred{i};
            end
            
        end

    end

end