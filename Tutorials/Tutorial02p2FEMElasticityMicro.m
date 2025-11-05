classdef Tutorial02p2FEMElasticityMicro < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
    end

    methods (Access = public)

        function obj = Tutorial02p2FEMElasticityMicro()
            obj.createMesh();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
        end

    end

    methods (Access = private)
        
        function createMesh(obj)
            fullmesh = UnitTriangleMesh(200,200);
            ls = obj.computeCircleLevelSet(fullmesh);
            sUm.backgroundMesh = fullmesh;
            sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);
            holeMesh = uMesh.createInnerMesh();
            obj.mesh = holeMesh;
        end

        function ls = computeCircleLevelSet(obj, mesh)
            % gPar.type          = 'Circle';
            % gPar.radius        = 0.25;
            % gPar.xCoorCenter   = 0.5;
            % gPar.yCoorCenter   = 0.5;
            % g                  = GeometricalFunction(gPar);
            % phiFun             = g.computeLevelSetFunction(mesh);
            % lsCircle           = phiFun.fValues;
            % ls = -lsCircle;
            % 
            % gPar.type           = 'FourPerpendicularBars';
            % gPar.barWidth       = 1;
            % gPar.leftBar_xMax   = 0.2;
            % gPar.rightBar_xMin  = 0.6;
            % gPar.bottomBar_yMax = 0.2;
            % gPar.topBar_yMin    = 0.6;

            gPar.type = 'DiagonalBars';
            gPar.leftBar_xMax = 0.35;   % right edge of left bar
            gPar.barWidth = 0.1;

            gPar.rightBar_xMin = 1 - gPar.leftBar_xMax;  % left edge of right bar
            gPar.bottomBar_yMax = gPar.leftBar_xMax ; % top edge of bottom bar
            gPar.topBar_yMin = gPar.rightBar_xMin;    % bottom edge of top bar            

            % gPar.type          = 'PerperndicularNFiber';
            % gPar.nFibers       = 5;
            % gPar.minxCoor      = 0;
            % gPar.maxxCoor      = 1;
            % gPar.minyCoor      = 0;
            % gPar.maxyCoor      = 1;
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            ls = lsCircle;
        end


      function computeElasticProperties(obj)
            E  = 1;
            nu = 1/3;
            obj.young   = ConstantFunction.create(E,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);
        end

        function createMaterial(obj)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = obj.young;
            s.poisson = obj.poisson;
            tensor    = Material.create(s);
            obj.material = tensor;
        end


        function solveElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MICRO';
            s.material = obj.material;
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            % Options: REDUCED-FLUC / MONOLITHIC-FLUC / MONOLITHIC-DISP
            s.solverCase = 'DIRECT';
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            fem.solve();
            obj.stateProblem = fem;
        end

        function bc = createBoundaryConditions(obj)
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
        end

    end

end