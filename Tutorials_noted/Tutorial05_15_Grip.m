classdef Tutorial05_15_Grip < handle

    properties (Access = private)
        filename 
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblemLoadBased
        physicalProblemMotionBased
        adjointProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        primalUpdater
        optimizer
        gJ

        J_MT
        J_MT_vector
        nGDI
        nMP
        Kp_bar
        Kp_bar_vector
        complianceFuncs
        motionBasedFuncs

        phi_imposed
        beta_adaptive
        J_current_monitor
        
    end

    methods (Access = public)

        function obj = Tutorial05_15_Grip()
               % obj.J_MT_vector = [-1, -2, -0.5]; % motion transmission. Input=1, J_MT corresponds to the output
                obj.Kp_bar_vector = [1 0.01 0.0001];
                

                for i = 1:length(obj.Kp_bar_vector)
                    obj.Kp_bar = obj.Kp_bar_vector(i);
                  %  obj.J_MT = obj.J_MT_vector(i);

                    fprintf('\n--- Starting optimization for K = %.4f (%d/%d) ---\n', ...
                    obj.Kp_bar, i, length(obj.Kp_bar_vector));
    
                    obj.nGDI = 1;
                    obj.nMP = 1;
                    %obj.Kp_bar = 0.01;
                 %   obj.phi_imposed = obj.J_MT;
                 %   obj.beta_adaptive = 0.01; 
                 %   obj.J_current_monitor = obj.J_MT;
                    
                    obj.init();
                    obj.createMesh();
                    obj.createDesignVariable();
                    obj.createFilter();
                    obj.createMaterialInterpolator();
    
                    obj.createElasticProblem();
                    obj.createComplianceFunctions();
                    obj.createMotionBasedFunctions();
    
                    obj.createVolumeConstraint();                
                    obj.createCost();
                    obj.createConstraint();
    
                    obj.createDualVariable();
                    obj.createPrimalUpdater();
                    obj.createOptimizer();
    
                    %obj.motionTransmissionAccuracy(); % Compute motion transmission accuracy
                    obj.printFinalDisplacement_v3(); % print document with the final FEM displacements
                    obj.printFinalDesignVariable(); % print final design variable
                    obj.saveFigures(); % save matlab figures (design variable and monitoring)

                    fprintf('--- Finished optimization for K = %.4f (%d/%d) ---\n', ...
                    obj.Kp_bar, i, length(obj.Kp_bar_vector));
                end
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj) % 2:1 mesh for the inverter
            % Generate coordinates
            x1 = linspace(-0.5,1.5,100);
            x2 = linspace(0,2.5,200);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'f');
            s.coord  = V(:,1:2);
            s.connec = F;
            %mesh = Mesh.create(s);
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) 0.5*ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s); 

            coords  = obj.mesh.coord;

            isLeft = coords(:,1) >= 0.3 & coords(:,1) <= 0.35 & coords(:,2) >= 1.35 & coords(:,2) <= 1.65  ;
            isRight = coords(:,1) >= 0.65 & coords(:,1) <= 0.70 & coords(:,2) >= 1.35 & coords(:,2) <= 1.65  ;
            isTop = coords(:,1) >= 0.35 & coords(:,1) <= 0.65 & coords(:,2) >= 1.65 & coords(:,2) <= 1.7;
            isBottom = coords(:,1) >= 0.35 & coords(:,1) <= 0.65 & coords(:,2) >= 1.30 & coords(:,2) <= 1.35;
            isGround = coords(:,2) <= 1e-8;
            isInput = coords(:,2) >= 0.6 & coords(:,2) <= 0.8 & coords(:,1) >= 0.4 & coords(:,1) <= 0.6;


            vals = dens.fun.fValues;
            vals(isInput | isLeft | isRight | isBottom | isTop | isGround) = 1;
            dens.fun.setFValues(vals);

            obj.designVariable = dens;

        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function f = createGradientFilter(obj)
            s.filterType = 'PDE';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1'); 
            f = Filter.create(s);
        end

        function createMaterialInterpolator(obj)
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
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function createElasticProblem(obj)
            % One FEM is needed per each GDI (input and output in this
            % case)
            obj.physicalProblemLoadBased = cell(obj.nGDI,1); % use cell to store different objects

            for i = 1:obj.nGDI % one iteration per active degree
                s.mesh = obj.mesh;
                s.scale = 'MACRO';
                s.material = obj.createMaterial();
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
                s.material = obj.createMaterial();
                s.dim = '2D';
                s.boundaryConditions = obj.createBoundaryConditionsMotionBased();
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
            s.volumeTarget = 0.5;
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
            gscale = 1; % Weight of the constraint, set to 1 in swan because the optcreateConstraintimizer has its own weight handling 

            for i = 1:obj.nMP
                cParams.MotionBasedStrainEnergy = obj.motionBasedFuncs{i};
                cParams.gscale = gscale;
                cParams.Kp_bar = obj.Kp_bar;
                s.shapeFunctions{i} = MotionBasedStiffnessConstraint(cParams);
            end

            s.shapeFunctions{obj.nMP + 1} = obj.volume;

            s.Msmooth = obj.createMassMatrix;
            obj.constraint = Constraint(s);
        end

        function createDualVariable(obj)
            nConstr = sum(obj.nMP)+1;
            s.nConstraints   = nConstr;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 500; % max step size, then etaNorm has to approve. If not approved, line search trial
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 500;
            s.tolerance      = 1e-8;
            nConstr = sum(obj.nMP)+1;
            s.constraintCase = repmat({'INEQUALITY'},1,nConstr);
            s.primalUpdater  = obj.primalUpdater;
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.1; % max design change per iter, default was 0.02
            s.gJFlowRatio    = 0.1; % "weight" of the constraint 0.1
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = ['InvDens'];
            s.applyNonDesignRegion = true;
            s.applySymmetry = false;
            %s.physicalProblem = obj.physicalProblem;

            coords   = obj.mesh.coord;

            isLeft = coords(:,1) >= 0.3 & coords(:,1) <= 0.35 & coords(:,2) >= 1.35 & coords(:,2) <= 1.65  ;
            isRight = coords(:,1) >= 0.65 & coords(:,1) <= 0.70 & coords(:,2) >= 1.35 & coords(:,2) <= 1.65  ;
            isTop = coords(:,1) >= 0.35 & coords(:,1) <= 0.65 & coords(:,2) >= 1.65 & coords(:,2) <= 1.7;
            isBottom = coords(:,1) >= 0.35 & coords(:,1) <= 0.65 & coords(:,2) >= 1.30 & coords(:,2) <= 1.35;
            isGround = coords(:,2) <= 1e-8;
            isInput = coords(:,2) >= 0.6 & coords(:,2) <= 0.8 & coords(:,1) >= 0.4 & coords(:,1) <= 0.6;

            s.nonDesignRegion = isInput | isLeft | isRight | isBottom | isTop | isGround;
            s.nonDesignValue  = 1;

            % s.iterCallBack =@(nIter) obj.adaptiveUpdate(nIter);
            % s.monitoringCallBack =@() [obj.phi_imposed; obj.J_current_monitor];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f{1},1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh = obj.mesh;
            m = Material.create(s);
        end

        function bc = createBoundaryConditionsLoadBased(obj, gdi_index) 
            isLeft =@(coor) coor(:,1) >= 0.3 & coor(:,1) <= 0.35 & coor(:,2) >= 1.35 & coor(:,2) <= 1.65  ;
            isRight =@(coor) coor(:,1) >= 0.65 & coor(:,1) <= 0.70 & coor(:,2) >= 1.35 & coor(:,2) <= 1.65  ;
            isTop =@(coor) coor(:,1) >= 0.35 & coor(:,1) <= 0.65 & coor(:,2) >= 1.65 & coor(:,2) <= 1.7;
            isBottom =@(coor) coor(:,1) >= 0.35 & coor(:,1) <= 0.65 & coor(:,2) >= 1.30 & coor(:,2) <= 1.35;
            isGround =@(coor) coor(:,2) <= 1e-8;
            isInput =@(coor) coor(:,2) >= 0.6 & coor(:,2) <= 0.8 & coor(:,1) >= 0.4 & coor(:,1) <= 0.6;

            % Side walls fixed (left fixed, midle only x fixed)
            sDir{1}.domain = isGround;
            sDir{1}.direction = [1,2];
            sDir{1}.value = 0;

            sPL{1}.domain = isInput;
            sPL{1}.direction = [2];
            sPL{1}.value = -1;
    
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
            isLeft =@(coor) coor(:,1) >= 0.3 & coor(:,1) <= 0.35 & coor(:,2) >= 1.35 & coor(:,2) <= 1.65  ;
            isRight =@(coor) coor(:,1) >= 0.65 & coor(:,1) <= 0.70 & coor(:,2) >= 1.35 & coor(:,2) <= 1.65  ;
            isTop =@(coor) coor(:,1) >= 0.35 & coor(:,1) <= 0.65 & coor(:,2) >= 1.65 & coor(:,2) <= 1.7;
            isBottom =@(coor) coor(:,1) >= 0.35 & coor(:,1) <= 0.65 & coor(:,2) >= 1.30 & coor(:,2) <= 1.35;
            isGround =@(coor) coor(:,2) <= 1e-8;
            isInput =@(coor) coor(:,2) >= 0.6 & coor(:,2) <= 0.8 & coor(:,1) >= 0.4 & coor(:,1) <= 0.6;

            sDir{1}.domain = isGround;
            sDir{1}.direction = [1,2];
            sDir{1}.value = 0;

            sDir{2}.domain = isInput;
            sDir{2}.direction = [2];
            sDir{2}.value = -0.2;

            sDir{3}.domain = isLeft;
            sDir{3}.direction = [1];
            sDir{3}.value = -0.1;

            sDir{4}.domain = isRight;
            sDir{4}.direction = [1];
            sDir{4}.value = 0.1;

            sDir{5}.domain = isTop;
            sDir{5}.direction = [2];
            sDir{5}.value = 0.1;

            sDir{6}.domain = isBottom;
            sDir{6}.direction = [2];
            sDir{6}.value = -0.1;


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
            titles = {'Compliance input GDI', 'Compliance output GDI'};

            for i = 1:obj.nGDI
                sC.mesh = obj.mesh;
                sC.stateProblem = obj.physicalProblemLoadBased{i};
                compFromTensor = ComplianceFromConstitutiveTensor(sC);

                s.title = titles{i};
                s.mesh = obj.mesh;
                s.filter = obj.filter;
                s.material = obj.createMaterial();
                s.complainceFromConstitutive = compFromTensor;
                obj.complianceFuncs{i} = ComplianceFunctional(s);
            end
        end
        
        function createMotionBasedFunctions(obj)
            % create a MotionBasedStrainEnergy per free MP
            obj.motionBasedFuncs = cell(obj.nMP,1);

            for i = 1:obj.nMP
                s.mesh = obj.mesh;
                s.filter = obj.filter;
                s.gradientFilter = obj.createGradientFilter();
                s.material = obj.createMaterial();
                s.stateProblem = obj.physicalProblemMotionBased{i};
                obj.motionBasedFuncs{i} = MotionBasedStrainEnergy(s);
            end
        end

        function motionTransmissionAccuracy(obj)
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
        

            isOutput = @(coor) coor(:,1) >= 0.9 & coor(:,1) <= 1.1 & coor(:,2) >= 1-1e-8;           
            outputNodes = isOutput(coor);
            u_out       = mean(uVals(outputNodes,2));
        
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

            % Unit displ of input
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

        function J_current = computeCurrentMT(obj)
            mat = obj.createMaterial();
            xD  = obj.designVariable.obtainDomainFunction();
            if iscell(xD); xD = xD{1}; end
            xF  = obj.filter.compute(xD, 2);
            mat.setDesignVariable({xF});
            C   = mat.obtainTensor();

            fem = obj.physicalProblemMotionBased{1};
            fem.updateMaterial(C);
            bc  = obj.createBoundaryConditionsMotionTransmission();
            fem.boundaryConditions = bc;
            fem.createBCApplier();
            fem.createSolver();
            fem.solve();

            uVals     = fem.uFun.fValues;
            coor      = obj.mesh.coord;
            isOutput  = @(c) c(:,1) >= 0.9 & c(:,1) <= 1.1 & c(:,2) >= 1-1e-8;
            outNodes  = isOutput(coor);
            u_out     = uVals(outNodes, 2);
            J_current = mean(u_out);

            if isnan(J_current)
                J_current = obj.J_MT;
            end
        end

        function adaptiveUpdate(obj,nIter)
            J_current = obj.computeCurrentMT();
            % if nIter < 50
            %     obj.phi_imposed = obj.J_MT;
            % else
                obj.phi_imposed = obj.phi_imposed + obj.beta_adaptive * (obj.J_MT - J_current);            obj.J_current_monitor = J_current;
            % end

            fem = obj.physicalProblemMotionBased{1};
            fem.boundaryConditions = obj.createBoundaryConditionsMotionBased();
            fem.createBCApplier();
            fem.createSolver();
        end


        function printFinalDisplacement_v2(obj) % works
            % --- EXTRACT YOUR FINAL DATA ---
            % Replace 'obj' with whatever your main framework object is called at the end
            nodes = obj.mesh.coord;     
            elements = obj.mesh.connec; 
            
            % Grab the raw vector from your saved variable
            global FINAL_DISPLACEMENT;
            u_raw = FINAL_DISPLACEMENT.fValues; 
            
            nNodes = size(nodes, 1);
            
            % --- RESHAPE THE DISPLACEMENT (BULLETPROOF STACKED) ---
            nNodes = size(nodes, 1);
            nComps = 2; % 2D problem
            expected_length = nNodes * nComps;
            
            u_raw_clean = u_raw(1 : expected_length);
            
            % The (:) forces them to stand up as vertical columns!
            u_X = u_raw_clean(1 : nNodes);
            u_X = u_X(:); 
            
            u_Y = u_raw_clean(nNodes + 1 : end);
            u_Y = u_Y(:);
            
            % Combine them side-by-side into a clean [N x 2] matrix
            u_phys = [u_X, u_Y];
            
            % Ensure nodes are 3D for ParaView
            if size(nodes, 2) == 2
                nodes = [nodes, zeros(nNodes, 1)];
            end
            
            % Ensure u_phys has a Z-component (add a column of zeros)
            if size(u_phys, 2) == 2
                u_phys = [u_phys, zeros(nNodes, 1)];
            end
            % ------------------------------------------------------
            
            % --- WRITE TO FILE ---
            fid = fopen('Final_Displacements.vtk', 'w');
            
            fprintf(fid, '# vtk DataFile Version 2.0\n');
            fprintf(fid, 'Custom FEM Results\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
            
            % Write Nodes
            fprintf(fid, 'POINTS %d float\n', nNodes);
            fprintf(fid, '%f %f %f\n', nodes');
            
            % Write Elements
            nElems = size(elements, 1);
            ptsPerElem = size(elements, 2);
            fprintf(fid, '\nCELLS %d %d\n', nElems, nElems * (ptsPerElem + 1));
            
            % Correct from MATLAB (1-based) to VTK (0-based) indexing
            elemData = [(ptsPerElem * ones(nElems, 1)), elements - 1];
            if ptsPerElem == 3      % Triangles
                fprintf(fid, '%d %d %d %d\n', elemData');
                cellType = 5; 
            elseif ptsPerElem == 4  % Quads
                fprintf(fid, '%d %d %d %d %d\n', elemData');
                cellType = 9; 
            end
            
            fprintf(fid, '\nCELL_TYPES %d\n', nElems);
            fprintf(fid, '%d\n', cellType * ones(nElems, 1));
            
            % Write Displacements
            fprintf(fid, '\nPOINT_DATA %d\n', nNodes);
            fprintf(fid, 'VECTORS Displacement float\n');
            fprintf(fid, '%f %f %f\n', u_phys');
            
            fclose(fid);
            disp('Successfully wrote Final_Displacements.vtk!');

        end
        
        function printFinalDisplacement_v3(obj)
            num_case = find(obj.Kp_bar_vector == obj.Kp_bar);
            namePrint = sprintf('D_grip_FinalDispl_J_%g',num_case);
            uFun = obj.physicalProblemMotionBased{1}.uFun;
            uFun.print(namePrint);
        end

        function saveFigures(obj)
            num_case = find(obj.Kp_bar_vector == obj.Kp_bar);

            fig_design = figure(1); 
            fig_monitor = figure(2);
           % fig_adaptive = figure(3);
            fig_monitor.WindowState = 'maximized';
          %  fig_adaptive.WindowState = 'maximized';
            drawnow;

            name_design = sprintf('D_Grip_DesignMap_Case_%g.png', num_case );
            name_monitor = sprintf('D_Grip_Monitoring_Case_%g.png', num_case);
            %name_adaptive = sprintf('D_Inv_Adaptive_Case_%g.png',   num_case);
            exportgraphics(fig_design, name_design, 'Resolution', 300);
            exportgraphics(fig_monitor, name_monitor, 'Resolution', 300);
          %  exportgraphics(fig_adaptive, name_adaptive, 'Resolution', 300);
            close all
        end

        function printFinalDesignVariable(obj)
            num_case = find(obj.Kp_bar_vector == obj.Kp_bar);
            namePrint = sprintf('D_Grip_DesignVariable_kCase_%g',num_case);
            obj.designVariable.fun.print(namePrint);
        end
     end
end