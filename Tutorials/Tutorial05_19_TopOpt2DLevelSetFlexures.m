classdef Tutorial05_19_TopOpt2DLevelSetFlexures < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        C
        dC
        physicalProblems
        cost
        constraint
        primalUpdater
        optimizer
        E0; nu0; E1; nu1

        doc % Selected/chosen DOCs, ["tx"]
        dof % Selected/chosen DOFs, ["rz"]
        emax % Max strain energy for DOFs
        mdoc % Vector with active DOC, [1 0 0]
        mdof % Vector with active DOF, [0 0 1]
        deg % Active degrees DOC+DOF, [1 0 1]
        ndeg % Number of active degrees, 2
        strainEnergyFuncs
    end

    methods (Access = public)

        function obj = Tutorial05_19_TopOpt2DLevelSetFlexures()
                % Degrees
                obj.doc = ["tx"]; %tx
                obj.dof = ["ty"]; % ty
                obj.emax = 0.5;

                obj.preprocessDegrees();
                
                obj.init();
                obj.createMesh();
                obj.createDesignVariable();
                obj.createFilter();
                obj.createMaterial();

                obj.createElasticProblem();
                obj.createStrainEnergyFunctions();

                obj.createCost();
                obj.createConstraint();

                obj.createPrimalUpdater();
                obj.createOptimizer();

                obj.printFinalDisplacement_v3(); % print document with the final FEM displacements
                obj.printFinalDesignVariable(); % print final design variable
                obj.saveFigures(); % save matlab figures (design variable and monitoring)
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.E0  = 1e-3;
            obj.nu0 = 1/3;
            obj.E1  = 1;
            obj.nu1 = 1/3;
        end

        function preprocessDegrees(obj)
            ldoc = length(obj.doc); % Number of DOCs
            ldof = length(obj.dof); % Number of DOFs

            % Initialize vectors of DOC and DOF
            obj.mdoc = zeros(ldoc,3);
            obj.mdof = zeros(ldof,3);
            
            degrees = ["tx","ty","rz"]; % possible DOC and DOF
            for i=1:ldoc; obj.mdoc(i,:) = strcmp(degrees,obj.doc(i)); end % Compare DOCs to degrees
            for i=1:ldof; obj.mdof(i,:) = strcmp(degrees,obj.dof(i)); end % Compare DOFs to degrees
            % Now, if doc="ty, tx", then mdoc=[0 1 0; 1 0 0]
            obj.mdoc = max(obj.mdoc,[],1) == 1; % For example, if doc="ty, tx", then mdoc=[1 1 0]
            obj.mdof = max(obj.mdof,[],1) == 1;
            
            obj.deg = find(max([obj.mdoc; obj.mdof])); % if mdoc=[0 1 0] and mdof=[1 0 0] then deg=[1 2] (indeces of active degrees)
            obj.ndeg = length(obj.deg); % Number of active degrees, used to avoid doing the FEA of an inactive degree
            
            % Checks that a degree is not both dof and doc, also that at
            % least 1
            % dof and 1 doc and no more than 2
            assert(all((obj.mdoc+obj.mdof) < 2),'overlap between DOC and DOF');
            assert(ldoc >= 1 & ldoc <= 2,'set of DOC too small/big');
            assert(ldof >= 1 & ldof <= 2,'set DOF too small/big');
            assert(obj.ndeg >= 2 & obj.ndeg <=3,'incorect no of rhs');
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(1,1,100,100);
        end

        function createDesignVariable(obj) % from tutorial 3
            s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterial(obj)
            N = obj.mesh.ndim;
            muA    = obj.computeMu(obj.E0,obj.nu0);
            kappaA = obj.computeKappa(obj.E0,obj.nu0,N);
            muB    = obj.computeMu(obj.E1,obj.nu1);
            kappaB = obj.computeKappa(obj.E1,obj.nu1,N);

            mu    = @(rho) SimpAllInterpolator.computeMu(muA,muB,kappaA,kappaB,rho,N); mu = @(rho) Expand(mu(rho),4);
            kappa  = @(rho) SimpAllInterpolator.computeKappa(muA,muB,kappaA,kappaB,rho,N); kappa = @(rho) Expand(kappa(rho),4);
            lambda = @(rho) kappa(rho) - (2/N)*mu(rho);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            obj.C  = @(rho) 2*mu(rho{1}).*I + lambda(rho{1}).*IxI;

            dmu     = @(rho) SimpAllInterpolator.computeMuDerivative(muA,muB,kappaA,kappaB,rho,N); dmu = @(rho) Expand(dmu(rho),4);
            dkappa  = @(rho) SimpAllInterpolator.computeKappaDerivative(muA,muB,kappaA,kappaB,rho,N); dkappa = @(rho) Expand(dkappa(rho),4);
            dlambda = @(rho) dkappa(rho) - (2/N)*dmu(rho);
            obj.dC  = @(rho) {2*dmu(rho{1}).*I + dlambda(rho{1}).*IxI};
        end

        function mu = computeMu(obj,E,nu)
            mu = E./(2*(1+nu));
        end

        function kappa = computeKappa(obj,E,nu,N)
            kappa = E./(N*(1-(N-1)*nu));
        end

        function createElasticProblem(obj)
            obj.physicalProblems = cell(obj.ndeg,1); % use cell to store different objects
            degrees_list = ["tx", "ty", "rz"];

            for i = 1:obj.ndeg % one iteration per active degree
                current_idx = obj.deg(i); % select the index of the current degree (from preprocess)
                current_degree = degrees_list(current_idx); % select the string degree

                s.mesh = obj.mesh;
                s.scale = 'MACRO';
                s.material = [];
                s.dim = '2D';
                s.boundaryConditions = obj.createBoundaryConditions(current_degree);
                s.interpolationType = 'LINEAR';
                s.solverType = 'REDUCED';
                s.solverMode = 'DISP'; % now we impose displacements
                s.solverCase = DirectSolver();
                fem = ElasticProblem(s);
                obj.physicalProblems{i} = fem; % store the FEM of this degree in the cell
            end

        end

        function c = createComplianceFromConstiutive(obj,i)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblems{i};
            c = ComplianceOnlyDirichlet(s);
        end

        function createStrainEnergyFunctions(obj)
            count = 1;
            for i = 1:obj.ndeg % One iteration per active degree
                deg_idx = obj.deg(i);
                if obj.mdoc(deg_idx) == 1 % If the degree is a doc
                    s.mesh                       = obj.mesh;
                    s.filter                     = obj.filter;
                    s.complainceFromConstitutive = obj.createComplianceFromConstiutive(i);
                    s.C                          = obj.C;
                    s.dC                         = obj.dC;
                    obj.strainEnergyFuncs{count} = ComplianceFunctional(s); % Store the strain energy object of the degree
                    count = count + 1;
                end
            end
        end

        function createCost(obj)
            % Initialize the functions to evaluate and their weights
            s.shapeFunctions = {};
            s.weights = [];

            ldoc = sum(obj.mdoc); % How many docs

            for i = 1:length(obj.strainEnergyFuncs) % One iteration per active deg
                s.shapeFunctions{i} = obj.strainEnergyFuncs{i}; % Take the FlexureStrainEnergy of the doc
                s.weights(i) = -1/ldoc; % Take the weight of the doc, with all docs having the same weight
            end

            s.Msmooth = obj.createMassMatrix();
            obj.cost  = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions = {};
            count = 1;

            for i = 1:obj.ndeg
                deg_idx = obj.deg(i);
                if obj.mdof(deg_idx) == 1 % Select the dofs
                    s.mesh                       = obj.mesh;
                    s.filter                     = obj.filter;
                    s.complainceFromConstitutive = obj.createComplianceFromConstiutive(i);
                    s.C                          = obj.C;
                    s.dC                         = obj.dC;
                    s.complianceTarget           = obj.emax;
                    s.shapeFunctions{count} = ComplianceConstraint(s);
                    count = count +1;
                end
            end
            s.Msmooth      = obj.createMassMatrix();
            obj.constraint = Constraint(s);
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
            s.maxIter        = 3;
            s.tolerance      = 1e-6; % default was 1e-8
            s.constraintCase = {'INEQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.delta          = 0.02; % max allowed change for level set(def 0.1)
            s.deltaMin       = 0.005; % default was 0.005
            s.etaStar       = 0.5; % weight for the constraints (def 0.2)
            s.etaMax0        = 1;
            s.etaMaxMin      = 0.01;
            s.gif            = false;
            s.gifName        = 'Tutorial05_14_LS';
            s.printing       = false;
            s.printName      = 'I_LS_';
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj, degree_type) 
            % Locate the walls of the domain
            xmin = min(obj.mesh.coord(:,1));
            xmax = max(obj.mesh.coord(:,1));
            ymin = min(obj.mesh.coord(:,2));
            ymax = max(obj.mesh.coord(:,2));

            % Length of the lower face
            Lx = xmax-xmin;

            % Select the upper and lower walls
            isBottom =@(coor) coor(:,2) <= ymin + 1e-8;
            isTop =@(coor) coor(:,2) >= ymax - 1e-8;

            % Bottom is fixed for tx, ty, rz
            sDir{1}.domain    = isBottom; % fixed
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            % Top nodes BC depend on the case
            if strcmp(degree_type, "tx")
                sDir{2}.domain    = isTop; % fixed
                sDir{2}.direction = [1];
                sDir{2}.value     = 1;

                sDir{3}.domain    = isTop; % fixed
                sDir{3}.direction = [2];
                sDir{3}.value     = 0;
            
            elseif strcmp(degree_type, "ty")
                sDir{2}.domain    = isTop; % fixed
                sDir{2}.direction = [1];
                sDir{2}.value     = 0;

                sDir{3}.domain    = isTop; % fixed
                sDir{3}.direction = [2];
                sDir{3}.value     = 1;

            elseif strcmp(degree_type, "rz")
                sDir{2}.domain    = isTop; % fixed
                sDir{2}.direction = [1];
                sDir{2}.value     = 1;

                top_logical = isTop(obj.mesh.coord); % Obtain which nodes are top wall
                top_coor    = obj.mesh.coord(top_logical, :); % coordinates of the top wall nodes

                sDir{3}.domain    = isTop; % fixed
                sDir{3}.direction = [2];
                sDir{3}.value     = 1-2*(top_coor(:,1)-xmin)/Lx; % expression for the rotation of the top wall
            end            


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
        
        function printFinalDisplacement_v3(obj)
            degrees_list = ["tx", "ty", "rz"];
            for i = 1:obj.ndeg
                deg_name = degrees_list(obj.deg(i));
                namePrint = sprintf('LS_Rot_FinalDispl_%s', deg_name);
                uFun = obj.physicalProblems{i}.uFun;
                uFun.print(namePrint);
            end
        end

        function saveFigures(obj)
          %  num_case = find(obj.k_vector == obj.k_case);
            fig_design = figure(1); 
            fig_monitor = figure(2);
            fig_monitor.WindowState = 'maximized';
            drawnow;
            name_design = sprintf('LS_Rot_DesVar.png' );
            name_monitor = sprintf('LS_Rot_Monitoring.png');
            exportgraphics(fig_design, name_design, 'Resolution', 300);
            exportgraphics(fig_monitor, name_monitor, 'Resolution', 300);
            close all
        end

        function printFinalDesignVariable(obj)
           % num_case = find(obj.k_vector == obj.k_case);
            namePrint = sprintf('LS_Rot_DesignVariable');
            obj.designVariable.fun.print(namePrint);
        end
     end
end