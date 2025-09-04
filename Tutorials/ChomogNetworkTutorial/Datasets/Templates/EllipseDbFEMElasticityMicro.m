classdef EllipseDbFEMElasticityMicro < handle
    properties (Access = public)
        mesh
        stateProblem
    end

    properties (Access = private)
        young
        poisson
        material
        gPar
    end

    methods (Access = public)

        function obj = EllipseDbFEMElasticityMicro(gPar)
            obj.gPar = gPar;
            obj.createMesh();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
        end

    end

    methods (Access = private)

        function createMesh(obj)
            fullmesh = UnitTriangleMesh(100,100);
            sUm.backgroundMesh = fullmesh;
            sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
            if obj.gPar.xSide ~= 0 && obj.gPar.ySide ~= 0

                uMesh = UnfittedMesh(sUm);
                ls = obj.computeEllipseLevelSet(fullmesh);
                uMesh.compute(ls);
                holeMesh = uMesh.createInnerMesh();
                obj.mesh = holeMesh;
            else
                obj.mesh = fullmesh;
            end
            
        end

        function ls = computeEllipseLevelSet(obj, mesh)
            obj.gPar.type        = 'SmoothRectangle';
            obj.gPar.xCoorCenter = 0.5;
            obj.gPar.yCoorCenter = 0.5;
            obj.gPar.pnorm       = 2.0;
            g                    = GeometricalFunction(obj.gPar);
            phiFun               = g.computeLevelSetFunction(mesh);
            lsEllipse            = phiFun.fValues;
            ls = -lsEllipse;
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