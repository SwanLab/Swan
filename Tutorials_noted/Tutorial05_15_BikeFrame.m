classdef Tutorial05_15_BikeFrame < handle

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
        gJFlowRatio_vector
        gJFlowRatio

        J_MT
        nGDI
        nMP
        Kp_bar
        Kp_bar_vector
        complianceFuncs
        motionBasedFuncs
        height % height of the domain (full)
        width % half of the domain (for symmetry)
        amplitude % amplitude of the sine/cosine MP
    end

    methods (Access = public)

        function obj = Tutorial05_15_BikeFrame()
                
                obj.nGDI = 1;
                obj.nMP = 0;
                %obj.Kp_bar_vector = [0.01 0.025]; 
                obj.gJFlowRatio = 0.3;

                obj.height = 1;
                obj.width = 2;
                %obj.amplitude = 0.05*obj.height;
                
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
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj) 
            % Generate coordinates
            x1 = linspace(0,obj.width,200);
            x2 = linspace(0,obj.height,100);
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
            s.fHandle = @(x) 0.1*ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s); 
            obj.designVariable = dens;

            coords   = obj.mesh.coord;
            % 
            % isSurface = coords(:,2) >= 0.97*obj.height;

            isTopLeft = coords(:,1) <= 0.03*obj.width & coords(:,2) >= obj.height-1e-8;
            isTopRight = coords(:,1) >= 0.75*obj.width & coords(:,1) <= 0.78*obj.width & coords(:,2) >= obj.height-1e-8;
            % 
            vals = dens.fun.fValues;
            vals(isTopLeft | isTopRight) = 1;
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
            %obj.physicalProblemMotionBased = cell(obj.nMP, 1);

            % for i = 1:obj.nMP
            %     s.mesh = obj.mesh;
            %     s.scale = 'MACRO';
            %     s.material = obj.createMaterial();
            %     s.dim = '2D';
            %     s.boundaryConditions = obj.createBoundaryConditionsMotionBased(i);
            %     s.interpolationType = 'LINEAR';
            %     s.solverType = 'REDUCED';
            %     s.solverMode = 'DISP'; 
            %     s.solverCase = DirectSolver();
            %     fem = ElasticProblem(s);
            %     obj.physicalProblemMotionBased{i} = fem; % store the FEM of this degree in the cell
            % end

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
            s.volumeTarget = 0.1;
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
            s.tauMax = 700; % max step size, then etaNorm has to approve. If not approved, line search trial
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 400;
            s.tolerance      = 1e-8;
            nConstr = sum(obj.nMP)+1;
            s.constraintCase = repmat({'INEQUALITY'},1,nConstr);
            s.primalUpdater  = obj.primalUpdater;
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.1; % max design change per iter, default was 0.02
            s.gJFlowRatio    = obj.gJFlowRatio; % "weight" of the constraint 0.1
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = ['InvDens'];
            s.applySymmetry = false;
            s.applyNonDesignRegion = true;
            %s.physicalProblem = obj.physicalProblem;

            % isSurface = obj.mesh.coord(:,2) >= 0.97*obj.height;
             isTopLeft = obj.mesh.coord(:,1) <= 0.03*obj.width & obj.mesh.coord(:,2) >= obj.height-1e-8;
            isTopRight = obj.mesh.coord(:,1) >= 0.75*obj.width & obj.mesh.coord(:,1) <= 0.78*obj.width & obj.mesh.coord(:,2) >= obj.height-1e-8;
            
            s.nonDesignRegion = isTopLeft | isTopRight;
            s.nonDesignValue  = 1;

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
            % The conditions correspond to half of the domain.
            isTopLeft =@(coor) coor(:,1) <= 0.03*obj.width & coor(:,2) >= obj.height-1e-8;
            isTopRight =@(coor) coor(:,1) >= 0.75*obj.width & coor(:,1) <= 0.78*obj.width & coor(:,2) >= obj.height-1e-8;
            isBottomLeft =@(coor) coor(:,1) >= 0.25*obj.width & coor(:,1) <= 0.28*obj.width & coor(:,2) <= 1e-8;
            isBottomRight =@(coor) coor(:,1) >= 0.97*obj.width & coor(:,2) <= 1e-8;

            % Side walls fixed
            % sDir{1}.domain = isTopLeft;
            % sDir{1}.direction = [1,2];
            % sDir{1}.value = 0;
            % 
            % sDir{2}.domain = isTopRight;
            % sDir{2}.direction = [1,2];
            % sDir{2}.value = 0;

            % sDir{1}.domain = isBottomRight;
            % sDir{1}.direction = [1,2];
            % sDir{1}.value = 0;

            sDir{1}.domain = isBottomLeft;
            sDir{1}.direction = [1,2];
            sDir{1}.value = 0;

            % Unit load at the GDI

            sPL{1}.domain = isTopRight;
            sPL{1}.direction = [2];
            sPL{1}.value = -1;

            sPL{2}.domain = isTopLeft;
            sPL{2}.direction = [2];
            sPL{2}.value = -1;

            sPL{3}.domain = isBottomRight;
            sPL{3}.direction = 2;
            sPL{3}.value = 1;

            sPL{4}.domain = isTopLeft;
            sPL{4}.direction = [1];
            sPL{4}.value = -0;

            sPL{5}.domain = isTopRight;
            sPL{5}.direction = [1];
            sPL{5}.value = 1;

            % sPL{3}.domain = isBottomRight;
            % sPL{3}.direction = [2];
            % sPL{3}.value = 1;
            % 
            % sPL{4}.domain = isBottomLeft;
            % sPL{4}.direction = [2];
            % sPL{4}.value = 2;



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
            titles = {'Compliance input S GDI', 'Compliance input C GDI','Compliance output GDI'};

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
            num_case_gJ = find(obj.gJFlowRatio_vector == obj.gJFlowRatio);
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