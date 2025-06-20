classdef ExploringOptimalShapeFromFusion < handle

    properties (Access = private)
        filename
        mesh       
        young
        poisson   
        material
        stateProblem
        compliance
        volume
        fractionVolume
        designVariable
        filter
        materialInterpolator
        
    end

    methods (Access = public)

        function obj = ExploringOptimalShapeFromFusion()
            obj.init()
            obj.createMesh();
            obj.createVolume();
            obj.createDesignVariable();
            obj.createFilter();
            obj.computeElasticProperties();
            %obj.createMaterialInterpolator();
            obj.createMaterial();
            obj.solveElasticProblem();
            %obj.createComplianceFromConstiutive();
            %obj.createCompliance();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'BEAM_3D_SF25';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end


        function createVolume(obj)
            obj.volume = obj.mesh.computeVolume();
            % Fa falta restar 
            boxVolume = 2e5;
            InitialVolume = 7.002e5;
            beamVolume = obj.volume-boxVolume;
            obj.fractionVolume = beamVolume/InitialVolume;
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            
            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end


        function computeElasticProperties(obj)
            E  = 71e3; %1; % canviar?
            nu = 1/3; % canviar?
            obj.young   = ConstantFunction.create(E,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);
        end

        function createMaterialInterpolator(obj)
            E0 = 10e3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 71e3;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end


        function createMaterial(obj)
            s.type                 = 'DensityBased';
            s.ptype                = 'ELASTIC';
            s.dim                  = '3D';%obj.mesh.ndim;
            s.dim                  = '3D';
            s.mesh                 = obj.mesh;
            s.density              = obj.designVariable;
            s.materialInterpolator = obj.materialInterpolator;
            s.young                = obj.young;
            s.poisson              = obj.poisson;
            tensor                 = Material.create(s);
            obj.material = tensor;
        end

        function solveElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material;
            s.dim = '3D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            fem.solve();
            obj.stateProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.stateProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.material; %obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
            [J,dJ] = c.computeFunctionAndGradient(obj.designVariable);
        end
        
        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            zMax    = max(obj.mesh.coord(:,3));
            isDir   = @(coor)  abs(coor(:,1))==20;
            isForce = @(coor)  (abs(coor(:,1))>0.995*xMax & abs(coor(:,2))>=0.25*yMax & abs(coor(:,2))<=0.75*yMax); % & abs(coor(:,3))>=0.749*zMax & abs(coor(:,3))<=0.75*zMax);

         
            s.fHandle = @(x)  (abs(x(1,:,:))>0.995*xMax & abs(x(2,:,:))>=0.25*yMax & abs(x(2,:,:))<=0.75*yMax); %  & abs(x(3,:,:))>=0.749*zMax & abs(x(3,:,:))<=0.75*zMax);
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            print(aFun.project('P1'),'Force 2','Paraview')
            
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 3;
            sPL{1}.value     = -10000;

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

      

    end


end