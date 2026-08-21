classdef Tutorial05_20_TopOpt2DDensity_MultiDofDeformableMirror < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        C
        dC
        physicalProblemLoadBased
        physicalProblemMotionBased
        volume
        cost
        constraint
        primalUpdater
        optimizer
        etaStarVector
        etaStar
        E0; nu0; E1; nu1

        J_MT
        nGDI
        nMP
        Kp_bar
        Kp_bar_vector
        complianceFuncs
        height % height of the domain (full)
        width % half of the domain (for symmetry)
        amplitude % amplitude of the sine/cosine MP
    end

    methods (Access = public)

        function obj = Tutorial05_20_TopOpt2DDensity_MultiDofDeformableMirror()
                
                obj.nGDI = 3;
                obj.nMP = 2;
                obj.Kp_bar_vector = [0.2 0.5]; 
                obj.etaStarVector = [0.7];

                obj.height = 1;
                obj.width = 1;
                obj.amplitude = 0.05*obj.height;
                
                for a=1:length(obj.etaStarVector)
                    obj.etaStar = obj.etaStarVector(a);
                    fprintf('\n--- Starting optimization for gJ = %.4f (%d/%d) ---\n', ...
                            obj.etaStar, a, length(obj.etaStarVector));
                    for i=1:length(obj.Kp_bar_vector)
                        obj.Kp_bar = obj.Kp_bar_vector(i);
                        fprintf('\n--- Starting optimization for Kp_bar = %.4f (%d/%d) ---\n', ...
                            obj.Kp_bar, i, length(obj.Kp_bar_vector));

                        obj.init();
                        obj.createMesh();
                        obj.createDesignVariable();
                        obj.createFilter();
                        obj.createMaterial();

                        obj.createElasticProblem();
                        obj.createComplianceFunctions();

                        obj.createVolumeConstraint();
                        obj.createCost();
                        obj.createConstraint();

                        obj.createPrimalUpdater();
                        obj.createOptimizer();

                        %obj.motionTransmissionAccuracy(); % Compute motion transmission accuracy
                        obj.printFinalDisplacement_v3(); % print document with the final FEM displacements
                        obj.printFinalDesignVariable(); % print final design variable
                        obj.saveFigures(); % save matlab figures (design variable and monitoring)

                        fprintf('--- Finished optimization for Kp_bar = %.4f (%d/%d) ---\n', ...
                            obj.Kp_bar, i, length(obj.Kp_bar_vector));
                    end
                    fprintf('\n--- Finished optimization for gJ = %.4f (%d/%d) ---\n', ...
                        obj.etaStar, a, length(obj.etaStarVector));
                end
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

        function createMesh(obj)
            obj.mesh = TriangleMesh(obj.width,obj.height,100,100);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = true;

            isOutput = @(x) x(:,2) >= 0.94*obj.height;
            isInputS = @(x) x(:,1) >= 0.95*obj.width & ...
                x(:,2) <= 0.25*obj.height & ...
                x(:,2) >= 0.15*obj.height;
            isInputC = @(x) x(:,1) >= 0.95*obj.width & ...
                x(:,2) <= 0.6*obj.height & ...
                x(:,2) >= 0.5*obj.height;
            s.isFixed = obj.computeFixedVolumeDomain(@(x) isOutput(x) | isInputS(x) | isInputC(x));

            dens = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function isFixed = computeFixedVolumeDomain(obj,cond)
            coor    = obj.mesh.coord;
            isFixed = find(cond(coor));
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
            % One FEM is needed per each GDI (input and output in this
            % case)
            obj.physicalProblemLoadBased = cell(obj.nGDI,1); % use cell to store different objects

            for i = 1:obj.nGDI % one iteration per active degree
                s.mesh = obj.mesh;
                s.scale = 'MACRO';
                s.material = [];
                s.dim = '2D';
                s.boundaryConditions = obj.createBoundaryConditionsLoadBased(i);
                s.interpolationType = 'LINEAR';
                s.solverType = 'REDUCED';
                s.solverMode = 'DISP'; 
                s.solverCase = DirectSolver();
                fem = ElasticProblem(s);
                obj.physicalProblemLoadBased{i} = fem; % store the FEM of this degree in the cell
            end

            % One fem per each free MP (only one for the inverter)
            obj.physicalProblemMotionBased = cell(obj.nMP, 1);

            for i = 1:obj.nMP
                s.mesh = obj.mesh;
                s.scale = 'MACRO';
                s.material = [];
                s.dim = '2D';
                s.boundaryConditions = obj.createBoundaryConditionsMotionBased(i);
                s.interpolationType = 'LINEAR';
                s.solverType = 'REDUCED';
                s.solverMode = 'DISP'; 
                s.solverCase = DirectSolver();
                fem = ElasticProblem(s);
                obj.physicalProblemMotionBased{i} = fem; % store the FEM of this degree in the cell
            end

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
            s.volumeTarget = 0.25;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end
        
        function createCost(obj)
            % Initialize the functions to evaluate and their weights
            s.shapeFunctions = {};
            s.weights = [];

            for i = 1:obj.nGDI
                s.shapeFunctions{i} = obj.complianceFuncs{i};
                s.weights(i) = 1;
            end

            s.Msmooth = obj.createMassMatrix;
            obj.cost = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions = {};

            for i = 1:obj.nMP
                cParams.mesh                       = obj.mesh;
                cParams.filter                     = obj.filter;
                cParams.complainceFromConstitutive = obj.createComplianceFromConstiutiveOnlyDirichlet(i);
                cParams.C                          = obj.C;
                cParams.dC                         = obj.dC;
                cParams.complianceTarget           = obj.Kp_bar;
                s.shapeFunctions{i} = ComplianceConstraint(cParams);
            end

            s.shapeFunctions{obj.nMP + 1} = obj.volume;

            s.Msmooth = obj.createMassMatrix;
            obj.constraint = Constraint(s);
        end

        function c = createComplianceFromConstiutiveOnlyDirichlet(obj,i)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblemMotionBased{i};
            c = ComplianceOnlyDirichlet(s);
        end

        function createPrimalUpdater(obj)
            rho      = obj.designVariable;
            fixedDof = rho.getFixedDofs();
            s.ub     = ones(size(rho.fun.fValues));
            s.lb     = zeros(size(rho.fun.fValues));
            s.lb(fixedDof) = 1;
            s.tauMax = 500;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 3;
            s.tolerance      = 1e-8;
            nConstr = sum(obj.nMP)+1;
            s.constraintCase = repmat({'INEQUALITY'},1,nConstr);
            s.primalUpdater  = obj.primalUpdater;
            s.ub             = 1;
            s.lb             = 0;
            s.delta          = 0.03; % max design change per iter, default was 0.02
            s.etaStar        = obj.etaStar; % "weight" of the constraint 0.1
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = ['InvDens'];
            s.applySymmetry = false;
            s.applyNonDesignRegion = true;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;

        end

        function bc = createBoundaryConditionsLoadBased(obj, gdi_index)
            % The conditions correspond to half of the domain.
            isWall =@(coor) coor(:,1) <= 1e-8 & coor(:,2) <= 0.5*obj.height;
            isMiddle =@(coor) coor(:,1) >= obj.width - 1e-8;
            
            mirror2inputLoadRelation = 2; % input load always = 1

            % Side walls fixed
            sDir{1}.domain = isWall;
            sDir{1}.direction = [1,2];
            sDir{1}.value = 0;

            sDir{2}.domain = isMiddle;
            sDir{2}.direction = [1];
            sDir{2}.value = 0;

            % Unit load at the GDI
            isInputS =@(coor) coor(:,1) >= 0.95*obj.width & coor(:,2) <= 0.25*obj.height & coor(:,2) >= 0.15*obj.height;
            isInputC =@(coor) coor(:,1) >= 0.95*obj.width & coor(:,2) <= 0.6*obj.height & coor(:,2) >= 0.5*obj.height;
            isOutput =@(coor) coor(:,2) >= 0.94*obj.height;

            isInputS_nodes = isInputS(obj.mesh.coord);
            nNodesInputS = sum(isInputS_nodes);

            isInputC_nodes = isInputC(obj.mesh.coord);
            nNodesInputC = sum(isInputC_nodes);

            isOutput_nodes = isOutput(obj.mesh.coord);
            nNodesOutput = sum(isOutput_nodes);

            if gdi_index == 1 % Input S
                sPL{1}.domain = isInputS;
                sPL{1}.direction = [2];
                sPL{1}.value = 1/nNodesInputS;
            elseif gdi_index == 2 % Input C
                sPL{1}.domain = isInputC;
                sPL{1}.direction = [2];
                sPL{1}.value = 1/nNodesInputC;
            elseif gdi_index == 3 % Output
                sPL{1}.domain = isOutput;
                sPL{1}.direction = [2];
                sPL{1}.value = -mirror2inputLoadRelation/nNodesOutput;
            else
                warning('Wrong GDI indeces');
            end        

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            
            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end

            s.dirichletFun = dirichletFun;
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    
        function bc = createBoundaryConditionsMotionBased(obj, mp_index) 
            % The conditions correspond to half of the domain.
            isWall =@(coor) coor(:,1) <= 1e-8 & coor(:,2) <= 0.5*obj.height;
            isMiddle =@(coor) coor(:,1) >= obj.width - 1e-8;

            % Side walls fixed
            sDir{1}.domain = isWall;
            sDir{1}.direction = [1,2];
            sDir{1}.value = 0;

            sDir{2}.domain = isMiddle;
            sDir{2}.direction = [1];
            sDir{2}.value = 0;

            isInputS =@(coor) coor(:,1) >= 0.95*obj.width & coor(:,2) <= 0.25*obj.height & coor(:,2) >= 0.15*obj.height;
            isInputC =@(coor) coor(:,1) >= 0.95*obj.width & coor(:,2) <= 0.6*obj.height & coor(:,2) >= 0.5*obj.height;
            isOutput =@(coor) coor(:,2) >= 0.94*obj.height;

            % Free MP corresponds to upwards input movement with downwards
            % output movement
            if mp_index == 1
                sDir{3}.domain    = isInputS;
                sDir{3}.direction = 2;
                sDir{3}.value     = 1;
                
                Output_logical = isOutput(obj.mesh.coord);
                Output_coor = obj.mesh.coord(Output_logical,:);
    
                sDir{4}.domain    = isOutput;
                sDir{4}.direction = 2;
                sDir{4}.value     = obj.amplitude * sin(Output_coor(:,1)/obj.width * 1.5*pi); % y(x)=A sin(x/w * 1.5*pi) 
            
                sDir{5}.domain = isInputC;
                sDir{5}.direction = [1,2];
                sDir{5}.value = 0;
            
            elseif mp_index == 2
                sDir{3}.domain    = isInputC;
                sDir{3}.direction = 2;
                sDir{3}.value     = 1;
                
                Output_logical = isOutput(obj.mesh.coord);
                Output_coor = obj.mesh.coord(Output_logical,:);
    
                sDir{4}.domain    = isOutput;
                sDir{4}.direction = 2;
                sDir{4}.value     = obj.amplitude * cos(Output_coor(:,1)/obj.width * pi); % y(x)=A cos(x/w * pi) 
          
                sDir{5}.domain = isInputS;
                sDir{5}.direction = [1,2];
                sDir{5}.value = 0;
            
            else
                warning('Wrong MP indeces');
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

        function createComplianceFunctions(obj)
            obj.complianceFuncs = cell(obj.nGDI,1);

            for i = 1:obj.nGDI
                sC.mesh = obj.mesh;
                sC.stateProblem = obj.physicalProblemLoadBased{i};
                compFromTensor = ComplianceFromConstitutiveTensor(sC);

                s.mesh = obj.mesh;
                s.filter = obj.filter;
                s.C = obj.C;
                s.dC = obj.dC;
                s.complainceFromConstitutive = compFromTensor;
                obj.complianceFuncs{i} = ComplianceFunctional(s);
            end
        end

        function motionTransmissionAccuracy(obj) % IDEALLY THE KP TARGET FOR COMPLIANCE IS ADAPTED WITH MT variable
            % Reuse existing motion-based FEM with updated BCs
            fem = obj.physicalProblemMotionBased{1};
            bc  = obj.createBoundaryConditionsMotionTransmission();
            fem.boundaryConditions = bc;
            fem.createBCApplier();
            fem.createSolver();
            fem.solve();
        
            % Read displacements
            uFun  = fem.uFun;
            coor  = obj.mesh.coord;
            uVals = uFun.fValues;
        
            % Select output nodes with density filter
            isOutput    = @(coor) coor(:,1) >= 0.9 & coor(:,1) <= 1.1 & coor(:,2) >= 1-1e-8;
            outputNodes = isOutput(coor);
            xD          = obj.designVariable.fun.fValues;
            outputDens  = xD(outputNodes);
            solidFilter   = outputDens > 0.5;
            u_out_all   = uVals(outputNodes, 2);
            u_out       = mean(u_out_all(solidFilter));
        
            % Motion transmission
            u_in = 1;
            MT   = u_out / u_in;
            eta  = 100 - abs(MT/obj.J_MT - 1) * 100;
        
            fprintf('\n--- Motion Transmission Accuracy ---\n');
            fprintf('J*    = %.4f (desired)\n', obj.J_MT);
            fprintf('J     = %.4f (obtained)\n', MT);
            fprintf('eta   = %.2f %%\n', eta);
            fprintf('------------------------------------\n');
        end

        function bc = createBoundaryConditionsMotionTransmission(obj)
            isLeft =@(coor) coor(:,1) <= 1e-8;
            isRight =@(coor) coor(:,1) >= 2-1e-8;

            % Side walls fixed
            sDir{1}.domain = isLeft;
            sDir{1}.direction = [1,2];
            sDir{1}.value = 0;

            sDir{2}.domain = isRight;
            sDir{2}.direction = [1,2];
            sDir{2}.value = 0;

            isInput =@(coor) coor(:,1) >= 0.9 & coor(:,1) <= 1.1 & coor(:,2) <= 1e-8;

            sDir{3}.domain = isInput;
            sDir{3}.direction = [2];
            sDir{3}.value = 1;    

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            
            pointloadFun = [];

            s.dirichletFun = dirichletFun;
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
        
        function printFinalDisplacement_v3(obj)
            num_case = find(obj.Kp_bar_vector == obj.Kp_bar);
            for i = 1:obj.nMP
                namePrint = sprintf('Mirror_FinalDispl_Kp_%g_MP_%g', num_case, i);
                uFun = obj.physicalProblemMotionBased{i}.uFun;
                uFun.print(namePrint);
            end

            for i=1:obj.nGDI
                namePrint = sprintf('Mirror_FinalDispl_Kp_%g_LP_%g', num_case, i);
                uFun = obj.physicalProblemLoadBased{i}.uFun;
                uFun.print(namePrint);
            end
        end

        function saveFigures(obj)
            num_case_gJ = find(obj.etaStarVector == obj.etaStar);
            num_case_Kp = find(obj.Kp_bar_vector == obj.Kp_bar);
            fig_design = figure(1); 
            fig_monitor = figure(2);
            fig_monitor.WindowState = 'maximized';
            drawnow;
            name_design = sprintf('Mirror_X_gJ_%g_Kp_%g.png', num_case_gJ, num_case_Kp );
            name_monitor = sprintf('Mirror_Monit_gJ_%g_Kp_%g.png', num_case_gJ, num_case_Kp);
            exportgraphics(fig_design, name_design, 'Resolution', 300);
            exportgraphics(fig_monitor, name_monitor, 'Resolution', 300);
            close all
        end

        function printFinalDesignVariable(obj)
            num_case = find(obj.Kp_bar_vector == obj.Kp_bar);
            namePrint = sprintf('D_Mirror_DesignVariable_kCase_%g',num_case);
            obj.designVariable.fun.print(namePrint);
        end
     end
end