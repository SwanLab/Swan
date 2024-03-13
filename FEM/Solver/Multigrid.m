classdef Multigrid < handle


    properties (Access = private)
        nDimf
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
        pdim
        RHS
        LHS
    end

    methods (Access = public)

        function obj = Multigrid(cParams)
            obj.init(cParams)
            s.type     = 'ELASTIC';
            s.scale    = 'MACRO';
            s.dim      = '2D';
            s.solverTyp = 'ITERATIVE';
            s.iterativeSolverType = 'CG';
            %obj.fem{1} = FEM.create(s);
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
            obj.LHS{1} = cParams.LHS;
            obj.RHS{1} = cParams.RHS;
            s.solverType = 'DIRECT';
            obj.solver{1} = Solver.create(s);
            obj.type     = cParams.type ;%'ELASTIC';
            obj.scale    = cParams.scale; %'MACRO';
            obj.pdim      = cParams.dim; %'2D';
            obj.nDimf    = cParams.nDimf;
        end

        function bc = createBoundaryConditions(obj,mesh,disp)
            rawBc    = obj.createRawBoundaryConditions(mesh);
            dim = getFunDims(obj,mesh,disp);
            rawBc.ndimf = dim.ndimf;
            rawBc.ndofs = dim.ndofs;
            s.mesh  = mesh;
            s.scale = 'MACRO';
            s.bc    = {rawBc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
%             obj.boundaryConditions{i+1} = bc;
        end
        
        function dim = getFunDims(obj,mesh,disp)
%             s.fValues = obj.mesh{i+1}.coord;
%             s.mesh = obj.mesh{i+1};
%             disp = P1Function(s);
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
            s.type = obj.type;
            s.scale = obj.scale;
            ngaus = 1;
            Id = ones(mesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = obj.pdim;
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

        function createFEMlevel(obj)
            for i = 1:obj.nLevel-1
%%%%%IN a new Class
                obj.createMatrixInterpolation(i)

               

                obj.mesh{i+1}     = obj.createMesh(i);
                obj.dispFun{i+1}  = P1Function.create(obj.mesh{i+1}, obj.nDimf);
                obj.boundaryConditions{i+1}       = obj.createBoundaryConditions(obj.mesh{i+1},obj.dispFun{i+1});
                obj.material{i+1} = obj.createMaterial(obj.mesh{i+1});
                s.solverTyp = 'ITERATIVE';
                s.iterativeSolverTyp = 'CG';
                s.tol                 = 1e-6;
                s.maxIter             = 20;

                obj.LHS{i+1} = obj.computeStiffnessMatrix(obj.mesh{i+1},obj.material{i+1},obj.dispFun{i+1});
                obj.RHS{i+1} = obj.createRHS(obj.mesh{i+1},obj.dispFun{i+1},obj.boundaryConditions{i+1});

            end
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

