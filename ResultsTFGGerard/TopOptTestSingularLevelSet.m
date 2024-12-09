classdef TopOptTestSingularLevelSet < handle

    properties (Access = private)
        filename %AFEGIDA PER LLEGIR MALLA
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TopOptTestSingularLevelSet()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'Malla_POCEXTENSA_6.m';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
            %%%%DENSITY^%%%%
            %s.fHandle = @(x) ones(size(x(1,:,:)));
            %s.ndimf   = 1;
            %s.mesh    = obj.mesh;
            %aFun      = AnalyticalFunction(s);
            %s.fun     = aFun.project('P1');
            %s.mesh    = obj.mesh;
            %s.type = 'Density';
            %s.plotting = false;
            %s.isFixed  = obj.computeFixedVolumeDomain(@(x) x(:,3)>=475, s.type);  %Volum no tocable
            %dens    = DesignVariable.create(s);
            %obj.designVariable = dens;
            %%%%LEVELSETv%%%%
             s.type = 'Full';
             g      = GeometricalFunction(s);
             lsFun  = g.computeLevelSetFunction(obj.mesh);
             s.fun  = lsFun;
             s.mesh = obj.mesh;
             s.type = 'LevelSet';
             s.plotting = false;
             s.isFixed  = obj.computeFixedVolumeDomain(@(x) x(:,3)>=475, s.type); % Jose: Aqui tambe es necessari
             ls     = DesignVariable.create(s);
             obj.designVariable = ls;
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

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialInterpolator(obj)
          %LEVEL-SET
            E0   = 1e-3;
            nu0  = 1/3;
            E1   = 1;
            nu1  = 1/3;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
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
             f = obj.filter.compute(f{1},1); 
             %Density: f = f.project('P1');
             s.type                 = 'DensityBased';
             s.density              = f;
             s.materialInterpolator = obj.materialInterpolator;
             s.dim                  = '3D';
             s.mesh                 = obj.mesh;
             m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
            %s.boundaryConditions = obj.createBoundaryConditionsGidAndMatlab(); %GiD per restriccions i Matlab per forces
            s.boundaryConditions = obj.createNewBoundaryConditionsWithGiD(); %obj.createBoundaryConditions();   %CANVI JOSE (new=GiD. altre=matlab)
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'CG';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstiutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.97;                               %VOLUM FINAL (volum target)
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
            nnodes  = obj.mesh.nnodes;
            indices = transpose(1:nnodes);
            vals    = ones(size(indices));                    %IGUAL PELS DOS (LEVELSET I DENSITY)
            h       = obj.mesh.computeMeanCellSize();
            M       = h^2*sparse(indices,indices,vals,nnodes,nnodes);   
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            %s.monitoring     = true;
            %s.cost           = obj.cost;
            %s.constraint     = obj.constraint;
            %s.designVariable = obj.designVariable;
            %s.dualVariable   = obj.dualVariable;
            %s.maxIter        = 2000;     %Iteracions
            %s.tolerance      = 1e-12;     %Hi havia 1e-8
            %s.constraintCase = {'EQUALITY'};
            %s.primal         = 'PROJECTED GRADIENT'; 
            %s.ub             = 1;
            %s.lb             = 0;
            %s.etaNorm        = 0.01; %HI HAVIA 0.001 (A menor valor menor oscilació del resultat?)
            %s.gJFlowRatio    = 0.35;   %hi havia 2   %major=complirconstraintrapid    menor=prioritzarminimitzarcost
            %opt = OptimizerNullSpace(s);
            %opt.solveProblem();
            %obj.optimizer = opt;
            %obj.designVariable.fun.print('Density_3D_Singular'); %Guarda la simulació automàticament per poder veure-la després a paraview
            %%%%%%%% Density  ^%%%%%%%%%

            %%%%%%%%%% LevelSet v %%%%%%%%
             s.monitoring     = true;
             s.cost           = obj.cost;
             s.constraint     = obj.constraint;
             s.designVariable = obj.designVariable;
             s.dualVariable   = obj.dualVariable;
             s.maxIter        = 5000;
             s.tolerance      = 1e-12;
             s.constraintCase = {'EQUALITY'};
             s.primal         = 'SLERP';
             s.etaNorm        = 0.01; %HI HAVIA 0.01 (Igual per això donava error?)
             s.gJFlowRatio    = 0.5;    %major=complirconstraintrapid    menor=prioritzarminimitzarcost
             opt = OptimizerNullSpace(s);
             opt.solveProblem();
             obj.optimizer = opt;
        end

        function bc = createBoundaryConditionsGidAndMatlab(obj)
%           xMax    = max(obj.mesh.coord(:,1));
%           yMax    = max(obj.mesh.coord(:,2));
%           zMax    = max(obj.mesh.coord(:,3));
            
          femReader = FemInputReader_GiD();
          s         = femReader.read(obj.filename);
          sDir      = obj.computeCondition(s.dirichlet);
          
          isForce1 = @(coor) coor(:,3)>=-225 & coor(:,3)<=240 & coor(:,2)>=47.2 & coor(:,2)<=47.4; %Eix Y comprovat amb paraview (valor exacte = 47.29999...)
          isForce2 = @(coor) coor(:,3)<=-225.1 & coor(:,3)>=-225.2; %Força aplicada a l'eix Z (valor exacte = 225.14999...)

          sPL{1}.domain    = @(coor) isForce1(coor);
          sPL{1}.direction = 2;
          sPL{1}.value     = -0.5;

          sPL{1}.domain    = @(coor) isForce2(coor);
          sPL{1}.direction = 3;
          sPL{1}.value     = -0.87; 

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

        function newbcGiD = createNewBoundaryConditionsWithGiD(obj)
            femReader = FemInputReader_GiD();
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
            newbcGiD = BoundaryConditions(s);
        end
    % 
    % end
    % 
    % methods (Static, Access=private)
        function sCond = computeCondition(obj,conditions)
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