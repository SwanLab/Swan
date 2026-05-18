classdef Tutorial05_15_TopOpt2DLevelSetFlexures < handle

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

        function obj = Tutorial05_15_TopOpt2DLevelSetFlexures()
                % Degrees
                obj.doc = ["tx"]; % tx, ty
                obj.dof = ["ty"]; % ty, rz
                obj.emax = 0.005; %1

                obj.preprocessDegrees();
                
                obj.init();
                obj.createMesh();
                obj.createDesignVariable();
                obj.createFilter();
                obj.createMaterialInterpolator();

                obj.createElasticProblem();
                obj.createStrainEnergyFunctions();

                obj.createVolumeConstraint();                
                obj.createCost();
                obj.createConstraint();

                obj.createDualVariable();
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
            
            % Checks that a degree is not both dof and doc, also that at least 1
            % dof and 1 doc and no more than 2
            assert(all((obj.mdoc+obj.mdof) < 2),'overlap between DOC and DOF');
            assert(ldoc >= 1 & ldoc <= 2,'set of DOC too small/big');
            assert(ldof >= 1 & ldof <= 2,'set DOF too small/big');
            assert(obj.ndeg >= 2 & obj.ndeg <=3,'incorect no of rhs');
        end

        function createMesh(obj) % must be modified
            % Generate coordinates
            x1 = linspace(0,1,100);
            x2 = linspace(0,1,100);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'f');
            s.coord  = V(:,1:2);
            s.connec = F;
            %mesh = Mesh.create(s);
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj) % from tutorial 3
            s.type = 'Holes';
            % For Holes
            s.dim = 2;
            s.nHoles = [80, 80];
            s.totalLengths = [1, 1];
            s.phases = [0, 0];
            s.phiZero = 0.3;
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
            obj.physicalProblem = cell(obj.ndeg,1); % use cell to store different objects
            degrees_list = ["tx", "ty", "rz"];

            for i = 1:obj.ndeg % one iteration per active degree
                current_idx = obj.deg(i); % select the index of the current degree (from preprocess)
                current_degree = degrees_list(current_idx); % select the string degree

                s.mesh = obj.mesh;
                s.scale = 'MACRO';
                s.material = obj.createMaterial();
                s.dim = '2D';
                s.boundaryConditions = obj.createBoundaryConditions(current_degree);
                s.interpolationType = 'LINEAR';
                s.solverType = 'REDUCED';
                s.solverMode = 'DISP'; % now we impose displacements
                s.solverCase = DirectSolver();
                fem = ElasticProblem(s);
                obj.physicalProblem{i} = fem; % stor the FEM of this degree in the cell
            end

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
            s.mesh           = obj.mesh
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

        function createStrainEnergyFunctions(obj)
            obj.strainEnergyFuncs = cell(obj.ndeg,1); % Initialize cell to store strain energy of each degree

            for i = 1:obj.ndeg % One iteration per active degree
                s.mesh         = obj.mesh;
                s.filter       = obj.filter;
                s.gradientFilter = obj.createGradientFilter;
                %s.material     = obj.materialInterpolator;
                s.material = obj.createMaterial();
                s.stateProblem = obj.physicalProblem{i};  % Give the FEM of that particular degree              
                obj.strainEnergyFuncs{i} = FlexureStrainEnergy(s); % Store the strain energy object of the degree
            end
        end
        

        function createCost(obj)
            % Initialize the functions to evaluate and their weights
            s.shapeFunctions = {};
            s.weights = [];

            ldoc = sum(obj.mdoc); % How many docs
            count = 1;

            for i = 1:obj.ndeg % One iteration per active deg
                deg_idx = obj.deg(i);
                if obj.mdoc(deg_idx) == 1 % If the degree is a doc
                    s.shapeFunctions{count} = obj.strainEnergyFuncs{i}; % Take the FlexureStrainEnergy of the doc
                    s.weights(count) = -1/ldoc; % Take the weight of the doc, with all docs having the same weight
                    count = count+1;
                end
            end

            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions = {};
            count = 1;
            gscale = 1; % Weight of the constraint, set to 1 in swan because the optimizer has its own weight handling 

            for i = 1:obj.ndeg
                deg_idx = obj.deg(i);
                if obj.mdof(deg_idx) == 1 % Select the dofs
                    cParams.strainEnergyFunc = obj.strainEnergyFuncs{i}; % take the strainEnergyFunc of that degree
                    cParams.gscale = gscale; % take its gscale
                    cParams.emax = obj.emax;
                    s.shapeFunctions{count} = FlexureDOFConstraint(cParams); % give the strainEnergyFunc and gscale to FlexureDOFConstraint to obtain the shape functions
                    count = count +1;
                end
            end

            %s.shapeFunctions{count} = obj.volume;

            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            nConstr = sum(obj.mdof);
            s.nConstraints   = nConstr;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createPrimalUpdater(obj) 
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj) % from mix tutorial 3
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 300;
            s.tolerance      = 1e-6; % default was 1e-8
            s.constraintCase = {'INEQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.1; % max allowed change for level set(def 0.1)
            s.etaNormMin     = 0.005; % default was 0.005
            s.gJFlowRatio    = 1; % weight for the constraints (def 0.2)
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
            s.gif            = false;
            s.gifName        = 'Tutorial05_14_LS';
            s.printing       = false;
            s.printName      = 'I_LS_';
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

        function bc = createBoundaryConditionsAdjoint(obj) % must be modified

            % the BC in the LS case have been modified, check if they are
            % the same for comparison purposes
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

            % NOTE: the order (index) of the input and output bc have been
            % switched so that they coincide with the indeces used when
            % computing the input and output costs separately

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
            namePrint = sprintf('D_Inv_FinalDispl_kCase_%g',num_case);
            uFun = obj.physicalProblem.uFun;
            uFun.print(namePrint);
        end

        function saveFigures(obj)
            num_case = find(obj.k_vector == obj.k_case);
            fig_design = figure(1); 
            fig_monitor = figure(2);
            fig_monitor.WindowState = 'maximized';
            drawnow;
            name_design = sprintf('D_Inv_DesignMap_kCase_%g.png', num_case );
            name_monitor = sprintf('D_Inv_Monitoring_kCase_%g.png', num_case);
            exportgraphics(fig_design, name_design, 'Resolution', 300);
            exportgraphics(fig_monitor, name_monitor, 'Resolution', 300);
            close all
        end

        function printFinalDesignVariable(obj)
            num_case = find(obj.k_vector == obj.k_case);
            namePrint = sprintf('D_Inv_DesignVariable_kCase_%g',num_case);
            obj.designVariable.fun.print(namePrint);
        end
     end
end