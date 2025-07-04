classdef ExploringOptimalShapeFromFusionCUBESAT < handle

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

        function obj = ExploringOptimalShapeFromFusionCUBESAT()
            obj.init()
            obj.createMesh();
            obj.createVolume();
            obj.createDesignVariable();
            obj.createFilter();
            obj.computeElasticProperties();
            obj.createMaterialInterpolator();
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
            file = 'CUBESAT_MALLA';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            m = s.mesh;
            con = m.connec;

            q = Quadrature.create(m,0);
            dv = m.computeDvolume(q);
            negElem=find(dv<=0);
            con(negElem,:) = [];
            sM.coord = m.coord;
            sM.connec = con;
            m2 = Mesh.create(sM);
            m2 = m2.computeCanonicalMesh();
            obj.mesh = m2;
        end


        function createVolume(obj)
            obj.volume = obj.mesh.computeVolume();
            % Fa falta restar 
            %boxVolume = 2e5;
            InitialVolume = 6.681e3;
           % preservedVolume = 
            beamVolume = obj.volume; %-boxVolume;
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
            E0 = 10;
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
            s.mesh                 = obj.mesh;
            s.density              = obj.designVariable;
            s.materialInterpolator = obj.materialInterpolator;
            % s.young                = obj.young;
            % s.poisson              = obj.poisson;
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
            s.solverCase = 'DIRECT'; % DIRECT
            fem = ElasticProblem(s);
            xD = obj.designVariable.obtainDomainFunction();
            obj.material.setDesignVariable(xD);
            C   = obj.material.obtainTensor();
            fem.updateMaterial(C);
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
            femReader = FemInputReaderGiD();
            s         = femReader.read(obj.filename);
            sPL       = obj.computeCondition(s.pointload);
            sDir      = obj.computeCondition(s.dirichlet);

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

    methods (Static, Access=private)
        function sCond = computeCondition(conditions)
            nodes = @(coor) 1:size(coor,1);
            dirs  = unique(conditions(:,2));
            j     = 0;
            for k = 1:length(dirs)
                rowsDirk = ismember(conditions(:,2),dirs(k));
                u        = unique(conditions(rowsDirk,3));
                for i = 1:length(u)
                    rows   = conditions(:,3)==u(i) & rowsDirk;
                    isCond = @(coor) ismember(nodes(coor),conditions(rows,1));
                    j      = j+1;
                    sCond{j}.domain    = @(coor) isCond(coor);
                    sCond{j}.direction = dirs(k);
                    sCond{j}.value     = u(i);
                end
            end
        end
    end


end