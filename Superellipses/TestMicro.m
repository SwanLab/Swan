classdef TestMicro < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
    end
    
    properties (Access = public)
        stateProblem
        geometryInclusion
    end

    methods (Access = public)

        function obj = TestMicro(cParams)
            obj.init(cParams);
            obj.createMesh();
            obj.mesh.plot();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
        end
        
        function init(obj,cParams)
            obj.geometryInclusion = cParams;
        end

    end

    methods (Access = private)
        
        function createMesh(obj)
            fullmesh = UnitQuadMesh(150,150);
            % fullmesh.plot();
            ls = obj.computeCircleLevelSet(fullmesh);
            sUm.backgroundMesh = fullmesh;
            sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);
            holeMesh = uMesh.createInnerMesh();
            obj.mesh = holeMesh;
        end

        function ls = computeCircleLevelSet(obj, mesh)
            gPar.type               = obj.geometryInclusion.type;
            % For the superellipse
            gPar.semiHorizontalAxis = obj.geometryInclusion.a;
            gPar.semiVerticalAxis   = obj.geometryInclusion.b;
            gPar.superEllipseFactor = obj.geometryInclusion.n;
            % For the superformula
%              gPar.semiHorizontalAxis = obj.geometryInclusion.a;
%              gPar.semiVerticalAxis   = obj.geometryInclusion.b;
%              gPar.m = obj.geometryInclusion.m;
%              gPar.n1 = obj.geometryInclusion.n1;
%              gPar.n2 = obj.geometryInclusion.n2;
%              gPar.n3 = obj.geometryInclusion.n3;
            gPar.xCoorCenter   = 0.5;
            gPar.yCoorCenter   = 0.5;
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            ls = -lsCircle;
        end


        function computeElasticProperties(obj)
            E1  = 1;
            nu1 = 1/3;
            E   = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu  = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            obj.young   = E;
            obj.poisson = nu;
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
