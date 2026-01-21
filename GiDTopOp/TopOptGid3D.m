classdef TopOptGid3D < handle

    properties (Access = private)
        filename
        mesh
        filter
        young
        poisson
        material
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
        primalUpdater
    end

    methods (Access = public)

        function obj = TopOptGid3D()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            %obj.computeElasticProperties();
            obj.createMaterialInterpolator();
            %obj.createMaterial();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();     
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)

            obj.mesh = HexaMesh(2,0.5,0.5,20,20,20);

            % file = 'BEAM_ANALYTIC';
            % obj.filename = file;
            % a.fileName = file;
            % s = FemDataContainer(a);
            % m = s.mesh;
            % con = m.connec;
            % 
            % q = Quadrature.create(m,0);
            % dv = m.computeDvolume(q);
            % negElem=find(dv<=0);
            % con(negElem,:) = [];
            % sM.coord = m.coord;
            % sM.connec = con;
            % m2 = Mesh.create(sM);
            % m2 = m2.computeCanonicalMesh();
            % obj.mesh = m2;
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            %bar1 = @(x) x(:,1)== 0;
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            %s.isFixed  = obj.computeFixedVolumeDomain(@(x) bar1(x), s.type);
            dens    = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
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
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f{1}.project('P1');            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '3D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
            
            % s.type                 = 'DensityBased';
            % s.ptype                = 'ELASTIC';
            % s.dim                  = '3D';%obj.mesh.ndim;
            % s.mesh                 = obj.mesh;
            % s.density              = obj.designVariable;
            % s.materialInterpolator = obj.materialInterpolator;
            % tensor                 = Material.create(s);
            % obj.material = tensor;
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'CG';
            fem = ElasticProblem(s);
            % xD = obj.designVariable.obtainDomainFunction();
            % obj.material.setDesignVariable(xD);
            % C   = obj.material.obtainTensor();
            % fem.updateMaterial(C);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
            %[J,dJ] = c.computeFunctionAndGradient(obj.designVariable);
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh,1,'P1');
            trial = LagrangianFunction.create(obj.mesh,1,'P1'); 
            M = IntegrateLHS(@(u,v) DP(v,u),test,trial,obj.mesh,'Domain');
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

       function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 100;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.01;
            s.gJFlowRatio    = 2;
            s.primalUpdater  = obj.primalUpdater;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            zMax    = max(obj.mesh.coord(:,3));
            isDir   = @(coor)  coor(:,1)==0;
            isForce = @(coor)  (coor(:,1)>0.995*xMax & abs(coor(:,2))>=0 & abs(coor(:,2))<yMax) & abs(coor(:,3))>=0 & abs(coor(:,3))<=zMax); % & abs(coor(:,3))>=0.749*zMax & abs(coor(:,3))<=0.75*zMax);
            
            numNodes = sum(isForce(obj.mesh.coord));
         
            s.fHandle = @(x)  (abs(x(1,:,:))>0.995*xMax & abs(x(2,:,:))>=0.25*yMax & abs(x(2,:,:))<=0.75*yMax); %  & abs(x(3,:,:))>=0.749*zMax & abs(x(3,:,:))<=0.75*zMax);
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

      %      print(aFun.project('P1'),'Force','Paraview')
            
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 3;
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end

        function isFixed = computeFixedVolumeDomain(obj,cond,type)
            coor  = obj.mesh.coord;
            nodes = find(cond(coor));
            isFixed.nodes = nodes;
            switch type
                case 'Density'
                    values = ones(size(nodes));
                case 'LevelSet'
                    values = -ones(size(nodes));
            end
            isFixed.values = values;
        end
    end
end