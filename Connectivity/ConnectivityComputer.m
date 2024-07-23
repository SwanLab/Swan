classdef ConnectivityComputer < handle

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        levelSet
        characteristicFunction
        designVariable

        materialInterpolation
        filter
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = ConnectivityComputer()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.createFilter();
            obj.createCharacteristicFunction();
            obj.createDesignVariable();

            % obj.levelSet.getUnfittedMesh().plot()
            % obj.density.plot()
            % shading flat
            % colormap('gray');
            % colormap(flipud(gray));
            % colorbar
            % figure

            obj.createCompliance();
            obj.computeEigenValueFunctional()
        end

    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function createMesh(obj)
            x1 = linspace(0,1,100);
            x2 = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createLevelSet(obj)
            s.type        = 'CircleInclusion';
            s.radius      = 0.4;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
            g             = GeometricalFunction(s);
            phi           = g.computeLevelSetFunction(obj.mesh);
            obj.levelSet = phi;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end        

        function createCharacteristicFunction(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh            = UnfittedMesh(s);
            uMesh.compute(obj.levelSet.fValues);

            obj.characteristicFunction  = CharacteristicFunction.create(uMesh);
        end

        function createDesignVariable(obj)
            s.fun  = obj.filter.compute(obj.characteristicFunction,3);
            s.mesh = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s);
            dens.plot();
            obj.designVariable = dens;
        end

        function createCompliance(obj)
            s.mesh                       = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            c.computeFunctionAndGradient(obj.designVariable);
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.createElasticProblem();
            c = ComplianceFromConstiutiveTensor(s);
        end

        function elas = createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';            
            s.interpolationType = 'LINEAR';
            elas = ElasticProblem(s);
        end

        function m = createMaterial(obj)
            dens = obj.designVariable.fun;
          %  f = x.obtainDomainFunction();
          %  f = f.project('P1');
            s.type                 = 'DensityBased';
            s.density              = dens;
            s.materialInterpolator = obj.createMaterialInterpolation();
            s.dim                  = '2D';
            m = Material.create(s);
        end        

        function m = createMaterialInterpolation(obj)
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
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end

        function computeEigenValueFunctional(obj)
            eigen = obj.computeEigenValueProblem();
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            mE = MinimumEigenValueFunctional(s);
            mE.computeFunctionAndGradient()
        end

        function eigen = computeEigenValueProblem(obj)
            s.mesh = obj.mesh;
            eigen  = StiffnessEigenModesComputer(s);
        end



    end

end
