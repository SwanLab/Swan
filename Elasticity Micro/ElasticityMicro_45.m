classdef ElasticityMicro_45 < handle

    properties (Access = private)
        young
        poisson
        material
        stateProblem
        mesh
        materialInterpolator
        filter
        designVariable
    end

    methods (Access = public)

        function obj = ElasticityMicro_45()
            obj.createMesh();
            obj.computeDensity();
            obj.computeElasticProperties();
            obj.createMaterialInterpolator(); 
            obj.solveElasticProblem();
            obj.computeVolume();
        end

    end

    methods (Access = private)
        
                function createMesh(obj)
            obj.mesh = UnitTriangleMesh(50,50);
        end

        function computeDensity(obj)
            [ls,phiFun] = obj.computeLevelSet(obj.mesh);
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);
            %holeMesh = uMesh.createInnerMesh();
            %obj.mesh = holeMesh;            
            close all;
            uMesh.plot;
            funLS     = CharacteristicFunction.create(uMesh);          
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s); 
            obj.filter = f;
            %obj.density = f.compute(funLS,2);
            s.fun = f.compute(funLS,2);
            s.type = 'Density';
            s.plotting = true;
            dens               = DesignVariable.create(s);
            obj.designVariable = dens;
            obj.designVariable.fun.plot
        end

        function [ls,phiFun] = computeLevelSet(obj, mesh)            
            g.type          = 'DiagonalNFibers';
            g.nFibers       = 4;
            g.minxCoor      = 0;
            g.maxxCoor      = 1;
            g.minyCoor      = 0;
            g.maxyCoor      = 1; 
            g                  = GeometricalFunction(g);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsValues           = phiFun.fValues;
            ls = lsValues;
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s); 
            obj.materialInterpolator = m;
        end

        function m = createMaterial(obj)          
            s.type                 = 'DensityBased';
            s.density              = obj.designVariable;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end        


      function computeElasticProperties(obj)
            E  = 1;
            nu = 1/3;
            obj.young   = ConstantFunction.create(E,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);
        end

        function solveElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MICRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            % Options: REDUCED-FLUC / MONOLITHIC-FLUC / MONOLITHIC-DISP
            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            s.density = obj.designVariable;
            s.filter = obj.filter;
            fem = ElasticProblemMicroAnisotropic(s);
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

        function computeVolume(obj)
           V  = Integrator.compute(obj.designVariable.fun,obj.mesh,2);
        end

    end

end