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
        filterConnect
        materialInterpolator
        minimumEigenValue
        lambdas
        phis
    end

    methods (Access = public)

        function obj = ExploringOptimalShapeFromFusion()
            obj.init()
            obj.createMesh();
            obj.createVolume();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createFilterConnectivity();
            obj.computeElasticProperties();
            obj.createMaterialInterpolator();
            obj.createMaterial();
            obj.solveElasticProblem();
            obj.computeEigenValue();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
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

        function createFilterConnectivity(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterConnect = f;
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

        function computeEigenValue(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect;
            %s.filterAdjoint     = obj.filterAdjointConnect;   
            s.targetEigenValue  = 50; % Minim eigenvalue      
            s.boundaryConditions = obj.createBoundaryConditions();
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
            n = 5;
            epsilon = 1e-5;
            p=8;
            [obj.lambda,obj.phis] = computeEigenValueFunctional(obj,n,epsilon,p);
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
            isDir   = @(coor)  coor(:,1)<=20*1.05;
            isForce = @(coor)  coor(:,1)>0.995*xMax;% & abs(coor(:,2))>=0.25*yMax & abs(coor(:,2))<=0.75*yMax); % & abs(coor(:,3))>=0.749*zMax & abs(coor(:,3))<=0.75*zMax);
            
            numNodes = sum(isForce(obj.mesh.coord));
         
            s.fHandle = @(x)  (abs(x(1,:,:))>0.995*xMax); % & abs(x(2,:,:))>=0.25*yMax & abs(x(2,:,:))<=0.75*yMax); %  & abs(x(3,:,:))>=0.749*zMax & abs(x(3,:,:))<=0.75*zMax);
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            print(aFun.project('P1'),'Force 2','Paraview')
            
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 3;
            sPL{1}.value     = -10000/numNodes;

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

         function  [lambdas, phis] = computeEigenValueFunctional(obj, n, epsilon, p)
            eigen = obj.computeEigenValueProblem(epsilon, p);
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            s.filter = obj.filterConnect;
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
%             mE = MinimumEigenValueFunctional(s);
            mE = MaximumEigenValueFunctional(s);
            [lambdas, phis] = mE.computeEigenModes(obj.designVariable, n);
         end

         function eigen = computeEigenValueProblem(obj,epsilon, p)
            s.mesh  = obj.mesh;
            s.epsilon = epsilon;
            s.p       = p;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            eigen   = StiffnessEigenModesComputer(s);
         end

         function  bc = createEigenvalueBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isFront = @(coor) abs(coor(:,2))==yMin;
            isBack = @(coor) abs(coor(:,2))== yMax;
            isDir   = @(coor) isLeft(coor) | isRight(coor) | isFront(coor) | isBack(coor);  
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = 0;
            sDir{1}.ndim = 1;
            
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);  
        end


      

    end


end