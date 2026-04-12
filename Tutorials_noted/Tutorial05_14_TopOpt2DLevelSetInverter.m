classdef Tutorial05_14_TopOpt2DLevelSetInverter < handle

    properties (Access = private)
        filename
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        adjointProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        primalUpdater
        optimizer
        gJ
        k_vector
        k_case
    end

    methods (Access = public)

        function obj = Tutorial05_14_TopOpt2DLevelSetInverter()
            obj.k_vector = [0.25 0.5 0.75];
            for a=1:length(obj.k_vector)               
                obj.k_case = obj.k_vector(a);
                fprintf('--- Starting Optimization for k = %f ---\n', obj.k_case);
                obj.init();
                obj.createMesh();
                obj.createDesignVariable();
                obj.createFilter();
                obj.createMaterialInterpolator();
                obj.createElasticProblem();
                obj.createAdjointProblem();
                obj.createNonSelfAdjCompliance();
                obj.createVolumeConstraint();
                obj.createCost();
                obj.createConstraint();
                obj.createDualVariable();
                obj.createPrimalUpdater();
                obj.createOptimizer();
                obj.printFinalDisplacement_v3(); % print document with the final FEM displacements
                obj.printFinalDesignVariable(); % print final design variable
                obj.saveFigures(); % save matlab figures (design variable and monitoring)
                fprintf('--- Optimization finished for k = %f ---\n', obj.k_case);
            end
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            % Generate coordinates
            x1 = linspace(0,1,100);
            x2 = linspace(0,1,100);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            %mesh = Mesh.create(s);
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj) % from tutorial 3
            s.type = 'Full';
            % For holes
            s.dim = 2;
            s.nHoles = [8, 8];
            s.totalLengths = [1, 1];
            s.phases = [0, 0];
            s.phiZero = 0.5;
            %
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
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialInterpolator(obj) % from tutorial 3
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
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function createAdjointProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditionsAdjoint();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.adjointProblem = fem;
        end

        function createNonSelfAdjCompliance(obj)
            s.mesh           = obj.mesh;
            s.filter         = obj.filter;
            s.material       = obj.createMaterial();
            s.stateProblem   = obj.physicalProblem;
            s.adjointProblem = obj.adjointProblem;
            c = NonSelfAdjointComplianceFunctional(s);
            obj.compliance = c;
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
            s.volumeTarget = 0.3; 
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
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
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

        function createPrimalUpdater(obj) % from tutorial 3
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj) % from mix tutorial 3
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'INEQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.01; % max allowed change for level set
            s.etaNormMin     = 0.005;
            s.gJFlowRatio    = 0.2; % weight for the constraints
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
            s.gif            = false;
            s.gifName        = 'Tutorial05_14_LS';
            s.printing       = false;
            s.printName      = 'G_LS_';
            s.k_case         = obj.k_case;
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
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function bc = createBoundaryConditions(obj)
            isDir   = @(coor)  coor(:,1)>=0 & coor(:,1)<=0.05 & coor(:,2)<=1e-8 | coor(:,1)<=1 & coor(:,1)>=0.95 & coor(:,2)<=1e-8; % bottom corners

            isPLTop      = @(coor)  (coor(:,1) >= 0.45 & coor(:,1) <= 0.55 & coor(:,2) == 1 ); % top part of the domain (output)
            isPLBottom   = @(coor)  (coor(:,1) >= 0.45 & coor(:,1) <= 0.55 & coor(:,2) == 0 );% bottom part of the domain (input)

            sDir{1}.domain    = @(coor) isDir(coor); % fixed
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isPLBottom(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = +10; % upward force on the input

            % sPL{2}.domain    = @(coor) isPLTop(coor);
            % sPL{2}.direction = 2;
            % sPL{2}.value     = 0;


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

        function bc = createBoundaryConditionsAdjoint(obj)
            
            isDir   = @(coor)  coor(:,1)>=0 & coor(:,1)<=0.05 & coor(:,2)<=1e-8 | coor(:,1)<=1 & coor(:,1)>=0.95 & coor(:,2)<=1e-8; % bottom corners

            isPLTop      = @(coor)  (coor(:,1) >= 0.45 & coor(:,1) <= 0.55 & coor(:,2) == 1 ); % top part of the domain (output)
            isPLBottom   = @(coor)  (coor(:,1) >= 0.45 & coor(:,1) <= 0.55 & coor(:,2) == 0 );% bottom part of the domain (input)

            sDir{1}.domain    = @(coor) isDir(coor); % fixed
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isPLBottom(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = +1; % upward force on the input

            sPL{2}.domain    = @(coor) isPLTop(coor);
            sPL{2}.direction = 2;
            sPL{2}.value     = +obj.k_case; % dummy load at the output



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

        function printFinalDisplacement(obj) % does not work
            % --- EXTRACT YOUR FINAL DATA ---
            % Replace 'obj' with whatever your main framework object is called at the end
            nodes = obj.mesh.coord;     
            elements = obj.mesh.connec; 
            
            % Grab the raw vector from your saved variable
            global FINAL_DISPLACEMENT;
            u_raw = FINAL_DISPLACEMENT.fValues; 
            
            nNodes = size(nodes, 1);
            
            % --- RESHAPE THE DISPLACEMENT (ROBUST VERSION) ---            
            nNodes = size(nodes, 1);
            nComps = 2; % It is a 2D tutorial, so we have X and Y components
            
            % 1. Calculate the exact number of displacement values we need
            expected_length = nNodes * nComps;
            
            % 2. Chop off any extra Lagrange multipliers at the end of the vector
            u_raw_clean = u_raw(1 : expected_length);
            
            % 3. Reshape safely into an [N x 2] matrix
            u_phys = reshape(u_raw_clean, nComps, nNodes)';
            
            % 4. Ensure nodes are 3D for ParaView
            if size(nodes, 2) == 2
                nodes = [nodes, zeros(nNodes, 1)];
            end
            
            % 5. Ensure u_phys has a Z-component (add a column of zeros)
            if size(u_phys, 2) == 2
                u_phys = [u_phys, zeros(nNodes, 1)];
            end
            % -------------------------------------------------
            
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
            num_case = find(obj.k_vector == obj.k_case);
            namePrint = sprintf('LS_FinalDispl_kCase_%g',num_case);
            uFun = obj.physicalProblem.uFun;
            uFun.print(namePrint);
        end

        function saveFigures(obj)
            num_case = find(obj.k_vector == obj.k_case);
            fig_design = figure(1); 
            fig_monitor = figure(2);
            fig_monitor.WindowState = 'maximized';
            drawnow;
            name_design = sprintf('LS_DesignMap_kCase_%g.png', num_case );
            name_monitor = sprintf('LS_Monitoring_kCase_%g.png', num_case);
            exportgraphics(fig_design, name_design, 'Resolution', 300);
            exportgraphics(fig_monitor, name_monitor, 'Resolution', 300);
            close all
        end

        function printFinalDesignVariable(obj)
            num_case = find(obj.k_vector == obj.k_case);
            namePrint = sprintf('LS_DesignVariable_kCase_%g',num_case);
            obj.designVariable.fun.print(namePrint);
        end


        
    end

end