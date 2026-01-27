classdef ParametrizedHomogenizationComputer < handle

    properties (Access = public)

    end

    properties (Access = private)
        backgroundMesh
        mesh
        material
        boundaryConditions
        stateProblem
    end

    methods (Access = public)

        function obj = ParametrizedHomogenizationComputer()
            obj.createBackgroundMesh();
            mx = linspace(0.02,0.95,20);
            my = linspace(0.02,0.95,20);
            for iMx = 1:length(mx)
                for iMy = 1:length(my)
                    obj.createMesh(mx(iMx),my(iMy));
                   % figure()
                   % obj.mesh.plot
                   % drawnow
                    obj.createMaterial();
                    obj.createBoundaryConditions();
                    Ch{iMx,iMy} = obj.solveElasticProblem();
                    close all
                  iMy
                end
                iMx
            end
            save('PaperDehomog/MicroStructureOptimization/OfflineRectangularChomog.mat',"Ch","mx","my")
        end
    end

    methods (Access = private)

        function createBackgroundMesh(obj)
            obj.backgroundMesh = UnitTriangleMesh(50,50);
            %obj.backgroundMesh = UnitQuadMesh(50,50);
        end
        
        function createMesh(obj,mx,my)
            ls = obj.computeCircleLevelSet(mx,my);
            s.backgroundMesh = obj.backgroundMesh;
            s.boundaryMesh   = obj.backgroundMesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(s);
            uMesh.compute(ls);
            holeMesh = uMesh.createInnerMesh();
            obj.mesh = holeMesh;
        end

        function lsV = computeCircleLevelSet(obj,mx,my)
            x = obj.backgroundMesh.coord(:,1);
            y = obj.backgroundMesh.coord(:,2);
            s.type    = 'RectangleInclusion';
            s.xSide = mx;
            s.ySide = my;
            s.xCoorCenter = (max(x) + min(x))/2;
            s.yCoorCenter = (max(y) + min(y))/2;
            g  = GeometricalFunction(s);
            ls = g.computeLevelSetFunction(obj.backgroundMesh);
            lsV = ls.fValues;
        end

        function createMaterial(obj)
            E  = 1;
            nu = 1/3;            
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = ConstantFunction.create(E,obj.mesh);
            s.poisson = ConstantFunction.create(nu,obj.mesh);
            tensor    = Material.create(s);
            obj.material = tensor;
        end


        function Ch = solveElasticProblem(obj)
            s.mesh     = obj.mesh;
            s.scale    = 'MICRO';
            s.material = obj.material;
            s.dim      = '2D';
            s.boundaryConditions = obj.boundaryConditions;
            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            fem.solve();
            Ch = fem.Chomog;
        end

        function createBoundaryConditions(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
            
            % Dirichlet
            
            sDir{1}.domain    = @(coor) isTop(coor) & isLeft(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
            
            sDir{2}.domain    = @(coor) isTop(coor) & isRight(coor);
            sDir{2}.direction = [1,2];
            sDir{2}.value     = 0;
            
            sDir{3}.domain    = @(coor) isBottom(coor) & isLeft(coor);
            sDir{3}.direction = [1,2];
            sDir{3}.value     = 0;
            
            sDir{4}.domain    = @(coor) isBottom(coor) & isRight(coor);
            sDir{4}.direction = [1,2];
            sDir{4}.value     = 0;
            
            % Periodic (NOTE: actually sets all boundaries as periodic hehe)
            
            sPer{1}.leader = @(coor) isLeft(coor);
            sPer{1}.follower = @(coor) isRight(coor);

            dirichletFun = [];
            periodicFun = [];

            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end

            for i = 1:numel(sPer)
                per = PeriodicCondition(obj.mesh, sPer{i});
                periodicFun = [periodicFun, per];
            end

            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = periodicFun;
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

    end

end