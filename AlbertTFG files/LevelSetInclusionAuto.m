classdef LevelSetInclusionAuto < handle
    
    properties (Access = public)
        stiffness
        strain
        stress
    end
    
    properties (Access = private)
        radius
        nodeDirection
        physicalProblem
        bc
        boundaryConditions
        problemSolver
        forces
        uFun
        strainFun
    end

    properties  (Access = protected)
        mesh
        material
        displacementFun
    end


    methods (Access = public)

        function [obj, u, L] = LevelSetInclusionAuto(r, i)
            obj.init(r, i)
            obj.createMesh();
            
            %% New ugly chunk of code warning
            [u, L] = obj.doElasticProblemHere();
            
            z.mesh      = obj.mesh;
            z.fValues   = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            z.order     = 'P1';
            uFeFun = LagrangianFunction(z);

            %obj.physicalProblem.solve();

            
            %% 
            
        end

    end

    methods (Access = private)

        function init(obj, r, i)
            % close all;
            % clc;
            obj.radius = r;
            obj.nodeDirection = i;
        end

        function createMesh(obj)
            bgMesh   = obj.createReferenceMesh();
            lvSet    = obj.createLevelSetFunction(bgMesh);
            uMesh    = obj.computeUnfittedMesh(bgMesh,lvSet);
            obj.mesh = uMesh.createInnerMesh();
        end

        function mesh = createReferenceMesh(~)
             %UnitMesh better
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            mesh = Mesh.create(s);
        end

        function levelSet = createLevelSetFunction(obj,bgMesh)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0.5;
            sLS.yCoorCenter = 0.5;
            sLS.radius      = obj.radius;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(bgMesh);
            levelSet        = lsFun.fValues;
        end

        function uMesh = computeUnfittedMesh(~,bgMesh,levelSet)
            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function s = createElasticProblem(obj)
            s.mesh     = obj.mesh;
            s.scale    = 'MACRO';
            obj.createMaterial();
            s.material = obj.material;
            s.dim      = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType  = 'LINEAR';
            s.solverType         = 'MONOLITHIC';
            s.solverMode         = 'DISP';
            s.solverCase         = 'DIRECT';
            % fem = ElasticProblem(s);
            % obj.physicalProblem = fem;
        end

        function createMaterial(obj)
            [young,poisson] = obj.computeElasticProperties(obj.mesh);
            s.type       = 'ISOTROPIC';
            s.ptype      = 'ELASTIC';
            s.ndim       = obj.mesh.ndim;
            s.young      = young;
            s.poisson    = poisson;
            obj.material = Material.create(s);
            
        end

        function [young,poisson] = computeElasticProperties(~,mesh)
            E  = 1;
            nu = 1/3;
            young   = ConstantFunction.create(E,mesh);
            poisson = ConstantFunction.create(nu,mesh);            
        end

        function bc = createBoundaryConditions(obj)

            v                    = zeros(8,1);
            v(obj.nodeDirection) = 1;
            nRes                 = [1 1 2 2 3 3 4 4]*v;
            assignMatrix         = [2 1 0 0 0 0 0 0
                                    0 0 2 1 0 0 0 0
                                    0 0 0 0 2 1 0 0
                                    0 0 0 0 0 0 2 1
                                    1 2 1 2 1 2 1 2];

            vSimp       = assignMatrix*v;
            dirs        = cell(5,1);
            [dirs{:,1}] = deal([1, 2]);
            [dirs{nRes}]  = deal(vSimp(nRes));
            [dirs{end}]   = deal(vSimp(end));

            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            tol     = 1e-10;

            corner1 = @(coor) abs(coor(:,1)-xMin) <= tol & abs(coor(:,2)-yMin)<= tol;
            corner2 = @(coor) abs(coor(:,1)-xMax) <= tol & abs(coor(:,2)-yMin)<= tol;
            corner3 = @(coor) abs(coor(:,1)-xMin) <= tol & abs(coor(:,2)-yMax)<= tol;
            corner4 = @(coor) abs(coor(:,1)-xMax) <= tol & abs(coor(:,2)-yMax)<= tol;

            cornerVec = {corner1; corner2; corner3; corner4};


            sDir{1}.domain    = @(coor) cornerVec{1}(coor);
            sDir{1}.direction = cell2mat(dirs(1));
            sDir{1}.value     = 0;

            
            sDir{2}.domain    = @(coor) cornerVec{2}(coor);
            sDir{2}.direction = cell2mat(dirs(2));
            sDir{2}.value     = 0;

            sDir{3}.domain    = @(coor) cornerVec{3}(coor);
            sDir{3}.direction = cell2mat(dirs(3));
            sDir{3}.value     = 0;

            sDir{4}.domain    = @(coor) cornerVec{4}(coor);
            sDir{4}.direction = cell2mat(dirs(4));
            sDir{4}.value     = 0;

            sDir{5}.domain    = @(coor) cornerVec{nRes}(coor);
            sDir{5}.direction = cell2mat(dirs(end));
            sDir{5}.value     = 1;

%             sPL{1}.domain    = @(coor) isForce(coor);
%             sPL{1}.direction = 2;
%             sPL{1}.value     = -1;
            sPL = {};

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = cat(2, dirichletFun, dir);
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = cat(2, pointloadFun, pl);
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

        function [u, L] = doElasticProblemHere(obj)
            cParams = obj.createElasticProblem();
            obj.createDisplacementFunHere();
            obj.createBCApplyerHere(cParams);
            obj.createSolverHere(cParams)
            obj.computeStiffnessMatrixHere();
            obj.computeForcesHere(cParams);
            [u, L]  = obj.computeDisplacementHere();
            obj.computeStrainHere();
            obj.computeStressHere();
        end

        function createDisplacementFunHere(obj)
            obj.displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function createBCApplyerHere(obj, cParams)
            s.mesh = obj.mesh;
            s.boundaryConditions = cParams.boundaryConditions;
            obj.bc = BCApplier(s);
        end

        function createSolverHere(obj, cParams)
            sS.type      = cParams.solverCase;
            solver       = Solver.create(sS);
            p.solverType = cParams.solverType;
            p.solverMode = cParams.solverMode;
            p.solver     = solver;

            p.boundaryConditions = cParams.boundaryConditions;
            p.BCApplier          = obj.bc;
            obj.problemSolver    = ProblemSolver(p);
        end

        function computeStiffnessMatrixHere(obj)
            ndimf      = obj.displacementFun.ndimf;
            m.type     = 'ElasticStiffnessMatrix';
            m.mesh     = obj.mesh;
            m.test     = LagrangianFunction.create(obj.mesh,ndimf, 'P1');
            m.trial    = obj.displacementFun;
            m.material = obj.material;

            m.quadratureOrder = 2;
            lhs               = LHSintegrator.create(m);
            obj.stiffness     = lhs.compute();
        end

        function computeForcesHere(obj, cParams)
            n.type         = 'Elastic';
            n.scale        = 'MACRO';
            n.dim          = obj.getFunDimsHere();
            n.BC           = cParams.boundaryConditions;
            n.mesh         = obj.mesh;
            n.material     = obj.material;
            n.globalConnec = obj.mesh.connec;

            RHSint = RHSintegrator.create(n);
            rhs    = RHSint.compute();
            % Perhaps move it inside RHSint?
            if strcmp(cParams.solverType,'REDUCED')
                R          = RHSint.computeReactions(obj.stiffness);
                obj.forces = rhs+R;
            else
                obj.forces = rhs;
            end
            
        end

        function dim = getFunDimsHere(obj)
            d.ndimf     = obj.displacementFun.ndimf;
            d.nnodes    = size(obj.displacementFun.fValues, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim         = d;
        end

        function [u, L] = computeDisplacementHere(obj)
            o.stiffness = obj.stiffness;
            o.forces    = obj.forces;
            [u,L]       = obj.problemSolver.solve(o);
            z.mesh      = obj.mesh;
            z.fValues   = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            z.order     = 'P1';
            obj.uFun    = LagrangianFunction(z);
            uSplit      = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.fValues = uSplit;
        end

        function computeStrainHere(obj)
            quad = Quadrature.create(obj.mesh, 2);
            xV   = quad.posgp;
            obj.strainFun  = SymGrad(obj.displacementFun);
%             strFun       = strFun.obtainVoigtFormat();
            obj.strain     = obj.strainFun.evaluate(xV);
        end

        function computeStressHere(obj)
            quad            = Quadrature.create(obj.mesh, 2);
            xV              = quad.posgp;
            stressFun       = DDP(obj.material, obj.strainFun);
            stressFun.ndimf = obj.strainFun.ndimf;
            obj.stress      = stressFun.evaluate(xV);

        end

    end
end