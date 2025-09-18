classdef TopOptCUBESATLevelSet < handle

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

        function obj = TopOptCUBESATLevelSet()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.computeElasticProperties();
            obj.createMaterialInterpolator();
            obj.createMaterial();
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
            file = 'L_SHAPE';
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

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            % isFixed for CUBESAT
            
            % bar1 = @(x) x(:,1)>=-113.15 & x(:,1)<=-106.65 & x(:,2) >= -50 & x(:,2)<= -43.5;
            % bar2 = @(x) x(:,1)>=-113.15 & x(:,1)<=-106.65 & x(:,2) >= 43.5 & x(:,2)<= 50;
            % bar3 = @(x) x(:,1)>=106.65 & x(:,1)<=113.15 & x(:,2) >= -50 & x(:,2)<= -43.5;
            % bar4 = @(x) x(:,1)>=106.65 & x(:,1)<=113.15 & x(:,2) >= 43.5 & x(:,2)<= 50;

             % isFixed for L-SHAPE?
           

            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            %s.isFixed  = obj.computeFixedVolumeDomain(@(x) bar1(x) | bar2(x) | bar3(x) | bar4(x), s.type);  %Volum no tocable
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
            tensor                 = Material.create(s);
            obj.material = tensor;
        end

        function createElasticProblem(obj)
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
            xD = obj.designVariable.obtainDomainFunction();
            obj.material.setDesignVariable(xD);
            C   = obj.material.obtainTensor();
            fem.updateMaterial(C);
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
            s.material                    = obj.material; %obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
            %[J,dJ] = c.computeFunctionAndGradient(obj.designVariable);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;                       %VOLUM FINAL
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
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSIntegrator.create(s);
            M = LHS.compute;     
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 100;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.02;
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 0.2;
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
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