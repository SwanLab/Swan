classdef TopOptTestTutorial3DDensityNullSpace < handle

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

        function obj = TopOptTestTutorial3DDensityNullSpace()
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
          %  obj.mesh = HexaMesh(2,1,1,40,20,20); %MALLA GENERADA PEL MATLAB

            %INTRODUIM COM GENERA LA MALLA EL GiD

            file = 'Malla_CARTESIANA';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;  %faltava ficar el obj
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = false;
            s.isFixed  = obj.computeFixedVolumeDomain(@(x) x(:,3)>=475, s.type);  %Volum no tocable
            dens    = DesignVariable.create(s);
            obj.designVariable = dens;
            %%%%DENSITY^%%%%
            %%%%LEVELSETv%%%%
            % s.type = 'Full';
            % g      = GeometricalFunction(s);
            % lsFun  = g.computeLevelSetFunction(obj.mesh);
            % s.fun  = lsFun;
            % s.mesh = obj.mesh;
            % s.type = 'LevelSet';
            % s.plotting = false;
            % ls     = DesignVariable.create(s);
            % obj.designVariable = ls;
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
            %%%%%%%%%%%%%%%%%DENSITY ^%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%LEVEL SET v%%%%%%%%%%%%%%%%%%
        %     E0   = 1e-3;
        %     nu0  = 1/3;
        %     E1   = 1;
        %     nu1  = 1/3;
        %     ndim = 2;
        % 
        %     matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
        %     matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);
        % 
        %     matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
        %     matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);
        % 
        %     s.typeOfMaterial = 'ISOTROPIC';
        %     s.interpolation  = 'SIMPALL';
        %     s.dim            = '3D';
        %     s.matA = matA;
        %     s.matB = matB;
        % 
        %     m = MaterialInterpolator.create(s);
        %     obj.materialInterpolator = m;
        
        
        end
         
         function m = createMaterial(obj)
             x = obj.designVariable;
             f = x.obtainDomainFunction();
             f = f.project('P1');
             s.type                 = 'DensityBased';
             s.density              = f;
             s.materialInterpolator = obj.materialInterpolator;
             s.dim                  = '3D';
             m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
            %s.boundaryConditions = obj.createBoundaryConditions(); 
            s.boundaryConditions = obj.createNewBoundaryConditionsWithGiD(); %obj.createBoundaryConditions();   %CANVI JOSE (new=GiD. altre=matlab)
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
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
            s.volumeTarget = 0.4;                               %VOLUM FINAL
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
             % s.test  = LagrangianFunction.create(obj.mesh,1,'P1'); %CANVI JOSE (NO TOCAR)
             % s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
             % s.mesh  = obj.mesh;
             % s.type  = 'MassMatrix';
             % LHS = LHSintegrator.create(s);
             % M = LHS.compute;

            nnodes  = obj.mesh.nnodes;
            indices = transpose(1:nnodes);
            vals    = ones(size(indices));
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
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 250;                       %Iteracions
            s.tolerance      = 1e-8;     %Hi havia 1e-8
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT'; 
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 1; %HI HAVIA 0.001 (A menor valor menor oscilació del resultat?)
            s.gJFlowRatio    = inf;   %hi havia 2   %major=complirconstraintrapid    menor=prioritzarminimitzarcost
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
            %%%%%%%% Density  ^%%%%%%%%%

            %%%%%%%%%% LevelSet v %%%%%%%%
            % s.monitoring     = true;
            % s.cost           = obj.cost;
            % s.constraint     = obj.constraint;
            % s.designVariable = obj.designVariable;
            % s.dualVariable   = obj.dualVariable;
            % s.maxIter        = 250;
            % s.tolerance      = 1e-8;
            % s.constraintCase = {'EQUALITY'};
            % s.primal         = 'SLERP';
            % s.etaNorm        = 0.05; %HI HAVIA 0.05 (A menor valor menor oscilació del resultat?)
            % s.gJFlowRatio    = 1;    %major=complirconstraintrapid    menor=prioritzarminimitzarcost
            % opt = OptimizerNullSpace(s);
            % opt.solveProblem();
            % obj.optimizer = opt;

        end

        function bc = createBoundaryConditions(obj)
            %---------------------------------------------------%
            %---------------BOUNDARY CONDITIONS GiD-------------%
            % femReader = FemInputReader_GiD();  %Llegim malla GiD
            % s = femReader.read(obj.filename);  %Llegim malla GiD
            % 
            % [u,v]=unique(s.pointload(:,3));  %Creem un vector [u,v] que es igual a totes les files de la tercera columna de s.pointload
            %                                  % u= valors diferents de les forces aplicades en diferents punts. v= nº de forces diferents
            % for i = 1:length(v)
            %     rows    = find(s.pointload(:,3)==u(i)); %Creem "rows" = numero de nodes que tenen força aplicada (1,2,...)
            %     isForce = @(coor) s.pointload(rows,1);  %isForce guarda tots els valors dels nodes amb força aplicada
            %     sPL{i}.domain    = @(coor) isForce(coor);
            %     sPL{i}.direction = s.pointload(v(i),2); %La direcció de la força serà la segona columna de s.pointload
            %     sPL{i}.value     = u(i);                %El valor de la força serà la u(v)
            % end
            % 
            % isDir   = @(coor)  unique(s.dirichlet(:,1));  %llegeix les constraints de desplaçaments del GiD
            % sDir{i}.domain    = @(coor) isDir(coor);
            % sDir{i}.direction = [1,2,3]; % [1,2] 2D  /  [1,2,3] 3D  %restricció en aquest punts.
            % sDir{i}.value     = 0;

            %---------------BOUNDARY CONDITIONS GiD-------------%
            %---------------------------------------------------%
            %---------------------------------------------------%
            %-------------BOUNDARY CONDITIONS MATLAB------------%

            %---------------3D CANTIELEVER CASE-------------%
          xMax    = max(obj.mesh.coord(:,1));
          yMax    = max(obj.mesh.coord(:,2));
          zMax    = max(obj.mesh.coord(:,3));

          isDir  = @(coor)  abs(coor(:,1))==0;
          isForce = @(coor) abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax & abs(coor(:,3))>=0.4*zMax & abs(coor(:,3))<=0.6*zMax;

          sDir{1}.domain    = @(coor) isDir(coor); %punt esquerre
          sDir{1}.direction = [1,2,3]; %restricció vertical, horitzontal i 3r eix
          sDir{1}.value     = 0;  %desplaçament =0

          sPL{1}.domain    = @(coor) isForce(coor);
          sPL{1}.direction = 2;
          sPL{1}.value     = -1;
%---------------3D CANTIELEVER CASE-------------%



            %-------------BOUNDARY CONDITIONS MATLAB------------%
            %---------------------------------------------------%
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