classdef MultigridTesting2 < handle
    
    
    properties (Access = private)
        nDimf
        boundaryConditionsFine
        fineMaterial
        quad
        basisFvalues
        basisVec
        eigenVec
        functionType
        nbasis
        Kmodal
        Mmodal
        D
        L
        Lt
        FineK
        Lchol
        FineKred
        coarseMesh
        CoarseK
        CoarseKred
        coarseMaterial
        boundaryConditionsCoarse
        I
        fineMeshCoord
        fineMeshConnec
        fineMesh
        FineDispFun
        CoarseDispFun
        FineFred
        CoarseFred
        data
    end
    
    methods (Access = public)
        
        function obj = MultigridTesting2()
            close all;
            addpath(genpath(fileparts(mfilename('fullpath'))))
            obj.init()
            obj.createCoarseMesh()
            obj.createMatrixInterpolation()
            obj.createFineMesh()
            obj.createBoundaryConditionsFine()
            obj.createBoundaryConditionsCoarse()
            obj.createFineMaterial()
            obj.createCoarseMaterial()
            obj.computeFineKred();
            obj.computeCoarseKred();
            obj.computeFineFred();
            obj.computeCoarseFred();
            %obj.solverMG();
            obj.createData();

        end

        function r = getdata(obj)
            r = obj.data;
        end

        function r = getBC(obj)
            r{1} = obj.boundaryConditionsCoarse;
            r{2} = obj.boundaryConditionsFine;
        end

        function r = getMesh(obj)
            r{1} = obj.coarseMesh;
            r{2} = obj.fineMesh;
        end
    end

    methods (Access = private)
        
        function init(obj)
            obj.nDimf = 2;
            obj.nbasis = 20;
            obj.functionType = 'P1';
        end

        function createCoarseMesh(obj)
            numero1 = 10;
            numero2 = 10;
            % Generate coordinates
            x1 = linspace(0,2,numero1);
            x2 = linspace(0,1,numero2);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            obj.coarseMesh = Mesh(s);
        end

        function createMatrixInterpolation(obj)
            p = obj.coarseMesh.coord;
            t = obj.coarseMesh.connec;

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
            obj.I = T;
        end

        function createFineMesh(obj)
            % x = obj.coarseMesh.coord(:,1);
            % y = obj.coarseMesh.coord(:,2);
            % obj.fineMeshCoord(:,1) = obj.I*x;
            % obj.fineMeshCoord(:,2) = obj.I*y;
            obj.fineMeshCoord = obj.I * obj.coarseMesh.coord;
            obj.fineMeshConnec = delaunayn(obj.fineMeshCoord);
            s.coord = obj.fineMeshCoord;
            s.connec = obj.fineMeshConnec;
            obj.fineMesh = Mesh(s);
        end

        function createBoundaryConditionsFine(obj)
            rawBc    = obj.createRawBoundaryConditionsFine();
            dim = getFunDimsFine(obj);
            rawBc.ndimf = dim.ndimf;
            rawBc.ndofs = dim.ndofs;
            s.mesh  = obj.fineMesh;
            s.scale = 'MACRO';
            s.bc    = {rawBc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditionsFine = bc;
        end

        function dim = getFunDimsFine(obj)
            s.fValues = obj.fineMesh.coord;
            s.mesh = obj.fineMesh;
            disp = P1Function(s);
            d.ndimf  = disp.ndimf;
            d.nnodes = size(disp.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.fineMesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function bc = createRawBoundaryConditionsFine(obj)
            dirichletNodes = abs(obj.fineMesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.fineMesh.coord(:,1));
            isInRight = abs(obj.fineMesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.fineMesh.coord(:,2)-1.5) < 0.1;
            %forceNodes = isInRight & isInMiddleEdge;
            forceNodes = isInRight;
            nodes = 1:obj.fineMesh.nnodes;
            bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            nodesdir=size(nodes(dirichletNodes),2);
            bcDir(1:nodesdir,end+1) = 1;
            bcDir(nodesdir+1:end,end) = 2;
            bcDir(:,end+1)=0;
            bc.dirichlet = bcDir;
            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = -1;
        end

        function createBoundaryConditionsCoarse(obj)
            rawBc    = obj.createRawBoundaryConditionsCoarse();
            dim = getFunDimsCoarse(obj);
            rawBc.ndimf = dim.ndimf;
            rawBc.ndofs = dim.ndofs;
            s.mesh  = obj.coarseMesh;
            s.scale = 'MACRO';
            s.bc    = {rawBc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditionsCoarse = bc;
        end

        function dim = getFunDimsCoarse(obj)
            s.fValues = obj.coarseMesh.coord;
            s.mesh = obj.coarseMesh;
            disp = P1Function(s);
            d.ndimf  = disp.ndimf;
            d.nnodes = size(disp.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.coarseMesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function bc = createRawBoundaryConditionsCoarse(obj)
            dirichletNodes = abs(obj.coarseMesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.coarseMesh.coord(:,1));
            isInRight = abs(obj.coarseMesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.coarseMesh.coord(:,2)-1.5) < 0.1;
            %forceNodes = isInRight & isInMiddleEdge;
            forceNodes = isInRight;
            nodes = 1:obj.coarseMesh.nnodes;
            bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            nodesdir=size(nodes(dirichletNodes),2);
            bcDir(1:nodesdir,end+1) = 1;
            bcDir(nodesdir+1:end,end) = 2;
            bcDir(:,end+1)=0;
            bc.dirichlet = bcDir;
            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = -1;
        end

        function createFineMaterial(obj)
            s.mesh = obj.fineMesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            ngaus = 1;
            Id = ones(obj.fineMesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.fineMesh.nelem;
            s.mesh  = obj.fineMesh;
            s.kappa = .9107*Id;
            s.mu    = .3446*Id;
            mat = Material.create(s);
            mat.compute(s);
            obj.fineMaterial = mat;
        end

        function createCoarseMaterial(obj)
            s.mesh = obj.coarseMesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            ngaus = 1;
            Id = ones(obj.coarseMesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.coarseMesh.nelem;
            s.mesh  = obj.coarseMesh;
            s.kappa = .9107*Id;
            s.mu    = .3446*Id;
            mat = Material.create(s);
            mat.compute(s);
            obj.coarseMaterial = mat;
        end

        function computeFineKred(obj)
            obj.FineDispFun = P1Function.create(obj.fineMesh, obj.nDimf);
            obj.FineK    = obj.computeStiffnessMatrix(obj.fineMesh,obj.fineMaterial,obj.FineDispFun);
            obj.FineKred = obj.boundaryConditionsFine.fullToReducedMatrix(obj.FineK);
        end

        function computeCoarseKred(obj)
            obj.CoarseDispFun = P1Function.create(obj.coarseMesh, obj.nDimf);
            obj.CoarseK    = obj.computeStiffnessMatrix(obj.coarseMesh,obj.coarseMaterial,obj.CoarseDispFun);
            obj.CoarseKred = obj.boundaryConditionsCoarse.fullToReducedMatrix(obj.CoarseK);
        end

        function k = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.fun      = displacementFun;
            % s.test      = displacementFun;
            % s.trial      = displacementFun;
            s.material = material;
            lhs = LHSintegrator.create(s);
            k   = lhs.compute();
        end 

        function computeFineFred(obj)
            RHS  = obj.createRHS(obj.fineMesh,obj.FineDispFun,obj.boundaryConditionsFine);
            Fext = RHS.compute();
            obj.FineFred = obj.boundaryConditionsFine.fullToReducedVector(Fext);
        end

        function computeCoarseFred(obj)
            RHS  = obj.createRHS(obj.coarseMesh,obj.CoarseDispFun,obj.boundaryConditionsCoarse);
            Fext = RHS.compute();
            obj.CoarseFred = obj.boundaryConditionsCoarse.fullToReducedVector(Fext);
        end
        
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
        end

        function solverMG(obj)
            hasNotConverged = true;
            while hasNotConverged

                hasNotConverged = norm(r) > tol;
            end

        end

        function vCycle(obj,A,B,xsol,Acoarse,Bcoarse)
            n = length(B);
            x = zeros(n,1);
            r = B - A * x;
            %             z = ModalTesting.applyPreconditioner(M,r);
            % z = obj.applyPrecoditionerMultigrid(r);
            %             z=r-z;
            p = r;
            rzold = r' * r;
            iter = 0;

            hasNotConverged = true;
            multigrid = true;

            while hasNotConverged
                
                Ap = A * p;
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                if iter == 0 || hasPartiallyConverged == true 
                    ri = r;
                end
                %                 z = ModalTesting.applyPreconditioner(M,r)
                rznew = r' * r;

                %hasNotConverged = sqrt(rsnew) > tol;
                hasNotConverged = norm(r) > tol;
                hasPartiallyConverged = norm(r)/norm(ri) < 0.5;

                if hasPartiallyConverged && multigrid
                    z = obj.applyPrecoditionerMultigrid(r,Acoarse,Bcoarse,A,B);
                    r = z;
                    multigrid = false;
                end

                p = r + (rznew / rzold) * p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
%                 err(iter)=norm(x-xsol);
%                 errAnorm(iter)=((x-xsol)')*A*(x-xsol);
            end

        end

        function createData(obj)
            obj.data(1).p = obj.coarseMesh.coord;
            obj.data(2).p = obj.fineMesh.coord;
            obj.data(1).t = obj.coarseMesh.connec;
            obj.data(2).t = obj.fineMesh.connec;
            %data(1).e = 
            %data(2).e = 
            obj.data(1).T = obj.I;
            obj.data(1).R = obj.I';
            obj.data(1).A = obj.CoarseKred;
            obj.data(2).A = obj.FineKred;
            obj.data(1).b = obj.CoarseFred;
            obj.data(2).b = obj.FineFred;
        end

    end

end