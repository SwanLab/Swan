classdef Tutorial02p2FEMElasticityMicro < handle
   
    properties (SetAccess = private, GetAccess = public)
        mesh
        young
        poisson
        bulk
        shear
        material
        stateProblem
    end

    methods (Access = public)

        function obj = Tutorial02p2FEMElasticityMicro(length,type,E,v)
            obj.createMesh(length,type);
            obj.computeElasticProperties(E,v);
            obj.createMaterial();
            obj.solveElasticProblem();
        end

    end

    methods (Access = private)
        
        function createMesh(obj,length,type)
            fullmesh = UnitTriangleMesh(300,300);
            ls = obj.computeCircleLevelSet(fullmesh,length,type);
            if type ~= "Full"
                sUm.backgroundMesh = fullmesh;
                sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
                uMesh              = UnfittedMesh(sUm);
                uMesh.compute(ls);
                holeMesh = uMesh.createInnerMesh();
                obj.mesh = holeMesh;
            else
                obj.mesh = fullmesh;
            end
        end

        function ls = computeCircleLevelSet(obj, mesh,length,type)
            if type == "Circle"
                gPar.type          = 'Circle';
                gPar.radius        = length;
                gPar.xCoorCenter   = 0.5;
                gPar.yCoorCenter   = 0.5;
                g                  = GeometricalFunction(gPar);
                phiFun             = g.computeLevelSetFunction(mesh);
                lsCircle           = phiFun.fValues;
                ls = -lsCircle;
            elseif type == "Square"
                gPar.type          = 'Square';
                gPar.length        = length;
                gPar.xCoorCenter   = 0.5;
                gPar.yCoorCenter   = 0.5;
                g                  = GeometricalFunction(gPar);
                phiFun             = g.computeLevelSetFunction(mesh);
                lsSquare           = phiFun.fValues;
                ls = -lsSquare;
            elseif type == "Full"
                gPar.type          = 'Full';
                g                  = GeometricalFunction(gPar);
                phiFun             = g.computeLevelSetFunction(mesh);
                lsFull             = phiFun.fValues;
                ls = -lsFull;
            end
        end


        function computeElasticProperties(obj,E1,nu1)
            E   = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu  = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            obj.young   = E;
            obj.poisson = nu;
        end

        function createMaterial(obj)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.bulk = obj.bulk;
            s.shear = obj.shear;
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