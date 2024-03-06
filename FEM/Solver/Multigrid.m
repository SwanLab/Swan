classdef Multigrid < handle


    properties (Access = private)
        ndimf
        data

        nLevel
        mesh
        I
        material
        boundaryConditions
        dispFun
        Kred
        Fred
        fem
        tol
        solver
        type
        scale  
        dim
    end

    methods (Access = public)

        function obj = Multigrid(cParams)
            obj.init(cParams)
            obj.createCoarseMesh(1)
            s.mesh = obj.mesh{1};
            s.bc   = obj.createBoundaryConditions(s.mesh);
            %             obj.createBoundaryConditions(0)
            %             obj.createMaterial(0)
            s.material = obj.createMaterial(s.mesh);
            s.type     = 'ELASTIC';
            s.scale    = 'MACRO';
            s.dim      = '2D';
            s.solverTyp = 'ITERATIVE';
            s.iterativeSolverType = 'CG';
            obj.fem{1} = FEM.create(s);
            obj.createFEMlevel();
            %             obj.computeKred(0)
            %             obj.computeFred(0)


            obj.createData();

        end

        function r = getdata(obj)
            r = obj.data;
        end

        function r = getBC(obj)
            r = obj.boundaryConditions;
        end

        function r = getMesh(obj)
            r = obj.mesh;
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
            obj.tol   = cParams.tol;
            obj.nLevel = cParams.nLevel;
            obj.mesh{1} = cParams.mesh;
            obj.boundaryConditions{1} = cParams.bc;
            obj.material{1} = cParams.material;
            s.solverType = 'DIRECT';
            obj.solver{1} = Solver.create(s);
            obj.type     = cParams.type ;%'ELASTIC';
            obj.scale    = cParams.scale; %'MACRO';
            obj.dim      = cParams.dim; %'2D';
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

        function meshFine = createMesh(obj,i)
            meshCoarse = obj.fem{i}.getMesh();
            meshCoord = obj.I{i} * meshCoarse.coord;
            meshConnec = delaunayn(meshCoord);
            s.coord = meshCoord;
            s.connec = meshConnec;
            %             obj.mesh{i+1} = Mesh(s);
            meshFine = Mesh(s);
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

        function createData(obj)

            for i = 1:obj.nMesh
                obj.data(i).p = obj.mesh{i}.coord;
                obj.data(i).t = obj.mesh{i}.connec;
                if i < obj.nMesh
                    obj.data(i).T = obj.I{i};
                    obj.data(i).R = obj.I{i}';
                end
                obj.data(i).A = obj.Kred{i};
                obj.data(i).b = obj.Fred{i};
            end

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

        %         function computeFred(obj,i)
        %             RHS  = obj.createRHS(obj.mesh{i+1},obj.dispFun{i+1},obj.boundaryConditions{i+1});
        %             Fext = RHS.compute();
        %             obj.Fred{i+1} = obj.boundaryConditions{i+1}.fullToReducedVector(Fext);
        %         end

        function RHS = createRHS(~,mesh,dispFun,boundaryConditions)
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

        function createFEMlevel(obj)
            for i = 1:obj.nMesh-1

                obj.createMatrixInterpolation(i)

                mesh     = obj.createMesh(i);
                dispFun  = P1Function.create(obj.mesh{i+1}, obj.nDimf);
                s.bc       = obj.createRawBoundaryConditions(s.mesh);
                s.material = obj.createMaterial(s.mesh);
                s.type     = 'ELASTIC';
                s.scale    = 'MACRO';
                s.dim      = '2D';
                s.solverTyp = 'ITERATIVE';
                s.iterativeSolverTyp = 'CG';
                s.tol                 = 1e-6;
                s.maxIter             = 20;

                LHS = obj.computeStiffnessMatrix(obj,s.mesh,s.material,displacementFun);
                s.type     = 'ElasticStiffnessMatrix';
                s.mesh     = mesh;
                s.fun      = displacementFun;
                % s.test      = displacementFun;
                % s.trial      = displacementFun;
                s.material = material;
                lhs = LHSintegrator.create(s);
                k   = lhs.compute();
            end

            obj.fem{i+1} = FEM.create(s);
        end
    end

    %          function createFEMlevel(obj)
    %             for i = 1:obj.nMesh-1
    %                 obj.createMatrixInterpolation(i)
    %                 s.mesh     = obj.createMesh(i);
    %                 s.bc       = obj.createRawBoundaryConditions(s.mesh);
    %                 s.material = obj.createMaterial(s.mesh);
    %                 s.type     = 'ELASTIC';
    %                 s.scale    = 'MACRO';
    %                 s.dim      = '2D';
    %                 s.solverTyp = 'ITERATIVE';
    %                 s.iterativeSolverTyp = 'CG';
    %                 s.tol                 = 1e-6;
    %                 s.maxIter             = 20;
    %
    %                 obj.fem{i+1} = FEM.create(s);
    %             end
    %         end

end

