classdef TutorialShells < handle


    properties (Access=private)
        mesh
        young
        area
        shear
        inertia
        density
        poisson, rotation
        uFun
        thetaFun
        wFun
        bcU,bcT,bcW
        lhs,RHSq
        solverType
        type, values, fun
        bcCase
        zLayer, interfaceIndex
        A_tensor, B_tensor, D_tensor, H_tensor
        materialProperties, stressCase
        material
        naturalFrequencies
        modeShapes
        massModal, dampingModal, stiffnessModal, FModal
        dampingRatio
        shearCorrectionFactor
    end

    methods (Access = public)
        %% TutorialShells
        function obj = TutorialShells()
            close all; clc;

            % 1. Initialization
            obj.createMesh('wingShape')    % 'unitTriangle' // 'wingShape'
            obj.createMaterial()
            obj.createSolutionField()
            obj.solverType = 'REDUCED';

            problemType    = 'FREE_VIBRATIONS';
            % Options: 'STATIC' / 'FREE_VIBRATIONS' / 'FORCED_VIBRATIONS'

            % 2. Boundary Conditions and Assembly
            obj.createBoundaryConditions();     % BOUNDARY CONDITIONS
            LHS = obj.createLHS();              % Stiffness matrix
            obj.lhs = LHS;

            % --- PRE-COMPUTE DOF LENGTHS AND INDICES ---
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            dofFU = obj.computeFreeDofs(obj.bcU);
            dofFT = obj.computeFreeDofs(obj.bcT);
            dofFW = obj.computeFreeDofs(obj.bcW);

            % 3. Resolution
            if strcmpi(problemType, 'FREE_VIBRATIONS')
                %% FREE VIBRATIONS
                fprintf('\n===== SOLVING FREE VIBRATION PROBLEM =====\n');

                MLHS = obj.createMassLHS();

                % Solve eigenvalue problem
                nModes = 10;
                minsol = 3;
                fprintf('Computing %d modes...\n', nModes);

                [modes, lambda] = eigs(LHS, MLHS, nModes, 'smallestabs');
                omega_squared = diag(lambda);

                % Check for negative or complex eigenvalues (indicates problem)
                if any(real(omega_squared) < 0)
                    warning('⚠ Negative eigenvalues detected! Check boundary conditions.');
                end
                if any(imag(omega_squared) ~= 0)
                    warning('⚠ Complex eigenvalues detected! Check matrices.');
                end

                % Take real part and ensure positive
                omega_squared = real(omega_squared);
                omega_squared(omega_squared < 0) = 0;  % Set negatives to zero
                omega_n = sqrt(omega_squared);  % rad/s
                freq_Hz = omega_n / (2*pi);     % Hz

                % Sort by frequency
                [freq_Hz, imodes] = sort(freq_Hz, 'ascend');
                omega_n = omega_n(imodes);
                modes = modes(:, imodes);

                % Display frequencies
                fprintf('\n===== NATURAL FREQUENCIES =====\n');
                for i = 1:nModes
                    fprintf('Mode %2d: %10.4f Hz  (%10.4f rad/s)\n', i, freq_Hz(i), omega_n(i));
                end
                fprintf('================================\n\n');

                % Save for later reference
                obj.naturalFrequencies = freq_Hz;
                obj.modeShapes = modes;

                % Process first N modes for plotting
                solutionsToProcess = min(minsol, nModes);
                timeInstants = [];  % Not used for FREE_VIBRATIONS

            elseif strcmpi(problemType, 'FORCED_VIBRATIONS')
                %% FORCED VIBRATIONS
                fprintf('\n===== SOLVING FORCED VIBRATION PROBLEM =====\n');

                % ========== FORCING PARAMETERS ==========
                f_force = 50;
                omega_force = 2*pi * f_force;               % Forcing frequency
                damping_ratio = obj.dampingRatio;           % Damping ratio

                % ========== INITIAL CONDITIONS ==========
                % All 0
                initialDisp = zeros(nU+nTheta+nW,1);
                initialVel = initialDisp;

                dynamicForceCase = 'STEP';
                % 'SINUSOIDAL' / 'STEP'

                MLHS = obj.createMassLHS();
                RHS = obj.createRHS();

                nModes = 10;
                fprintf('Computing %d modes for modal superposition...\n', nModes);

                [modal_shapes, lambda] = eigs(LHS, MLHS, nModes, 'smallestabs');

                % Process eigenvalues
                omega_squared = (diag(lambda));
                if any(real(omega_squared) < 0)
                    warning('⚠ Negative eigenvalues detected! Check boundary conditions.');
                end
                if any(imag(omega_squared) ~= 0)
                    warning('⚠ Complex eigenvalues detected! Check matrices.');
                end
                omega_squared = real(omega_squared);
                omega_squared(omega_squared < 0) = 0;
                omega_n = sqrt(omega_squared);
                freq_Hz = omega_n / (2*pi);

                % Sort
                [freq_Hz, imodes] = sort(freq_Hz, 'ascend');
                omega_n = omega_n(imodes);
                modal_shapes = modal_shapes(:, imodes);

                % Set Newmark Method
                fmax = max([max(freq_Hz),f_force]);
                fsampling = 10 * fmax;
                dt = 1 / fsampling;
                tfinal = 10 / freq_Hz(1);
                time = 0:dt:tfinal;
                nt = length(time);

                fprintf('\n===== TIME INTEGRATION SETUP =====\n');
                fprintf('Time step (dt):     %.6f s\n', dt);
                fprintf('Final time:         %.4f s\n', tfinal);
                fprintf('Number of steps:    %d\n', nt);
                fprintf('==================================\n\n');

                dynamicDisp = zeros(nU+nTheta+nW,nt);
                dynamicVel = dynamicDisp;
                dynamicDisp(:,1) = initialDisp;
                dynamicVel(:,1) = initialVel;

                modalDisp = modal_shapes' * dynamicDisp;
                modalVel = modal_shapes' * dynamicVel;

                % Dynamic force
                switch dynamicForceCase
                    case 'SINUSOIDAL'
                        fdynamic = sin(omega_force*time);
                    case 'STEP'
                        t0 = 0;
                        fdynamic = heaviside(time-t0);
                    otherwise
                        fdynamic = ones(1, nt);
                end

                % Modal Force, Mass, Stiffness, Damping
                dynamicFModal = modal_shapes' * RHS;
                obj.FModal = repmat(dynamicFModal,1,nt);
                obj.FModal = obj.FModal .* fdynamic;
                obj.massModal = modal_shapes' * MLHS * modal_shapes;
                obj.stiffnessModal = modal_shapes' * LHS * modal_shapes;
                obj.dampingModal = damping_ratio * eye(size(obj.massModal));

                % Modal accelerations
                modalAcc = zeros(size(modalDisp));
                modalAcc(:,1) = obj.massModal \ (obj.FModal(:,1) - obj.stiffnessModal*modalDisp(:,1) - obj.dampingModal*modalVel(:,1));

                % Check for near-resonance conditions
                for i = 1:nModes
                    freq_diff = abs(omega_force/(2*pi) - freq_Hz(i)) / freq_Hz(i);
                    if freq_diff < 0.1  % Within 10%
                        warning('⚠ WARNING: Forcing frequency near mode %d resonance (%.2f Hz)\n', i, freq_Hz(i));
                    end
                end

                % Newmark Method
                fprintf('Starting Newmark-β time integration...\n');
                tic;
                [next_modalDisp, next_modalVel, next_modalAcc] = obj.NewmarkMethod(dt);

                for i = 1:nt-1
                    modalDisp(:,i+1) = next_modalDisp(modalDisp(:,i), modalVel(:,i), modalAcc(:,i), obj.FModal(:,i+1));
                    modalAcc(:,i+1) = next_modalAcc(modalDisp(:,i+1), modalDisp(:,i), modalVel(:,i), modalAcc(:,i));
                    modalVel(:,i+1) = next_modalVel(modalVel(:,i), modalAcc(:,i), modalAcc(:,i+1));
                end
                elapsed = toc;
                fprintf('✓ Time integration complete in %.2f seconds\n', elapsed);

                % Reconstruct physical displacements for all time steps
                modes = modal_shapes * modalDisp;  % (nDOF x nt)

                % ========== SELECT TIME INSTANTS FOR POST-PROCESSING ==========
                % Option 1: Process all time steps (WARNING: can be slow!)
                % timeInstants = 1:nt;

                % Option 2: Process specific time instants
                % Examples:
                % - Final time only
                % timeInstants = nt;

                % - Every N-th timestep
                % N_skip = 10;
                % timeInstants = 1:N_skip:nt;

                % - Specific times (e.g., peak response, steady state)
                % Find max displacement

                w_history = modes(nU+nTheta+1:end, :);
                [max_w_value, linear_idx] = max(abs(w_history(:)));
                [node_idx_local, time_idx] = ind2sub(size(w_history), linear_idx);



                % timeInstants = [1, round(nt/4), round(nt/2), round(3*nt/4), idx_max, nt];
                % timeInstants = unique(timeInstants);  % Remove duplicates

                % timeInstants = [round(nt/4), idx_max, nt];
                timeInstants = [time_idx];
                timeInstants = unique(timeInstants);

                solutionsToProcess = length(timeInstants);

                fprintf('\n===== POST-PROCESSING SETUP =====\n');
                fprintf('Total time steps computed: %d\n', nt);
                fprintf('Time instants to process:  %d\n', solutionsToProcess);
                fprintf('Times: ');
                for i = 1:solutionsToProcess
                    fprintf('%.4f ', time(timeInstants(i)));
                end
                fprintf('s\n');
                fprintf('=================================\n\n');

                obj.animateDynamicResponse(modes, time, nU, nTheta, nW)


            elseif strcmpi(problemType, 'STATIC')
                %% STATIC
                fprintf('\n===== SOLVING STATIC PROBLEM =====\n');
                RHS = obj.createRHS();
                x = LHS \ RHS;
                modes = x;
                solutionsToProcess = 1;
                timeInstants = [];  % Not used for STATIC

            else
                error('Unknown problem type: %s. Use ''STATIC'', ''FREE_VIBRATIONS'', or ''FORCED_VIBRATIONS''.', problemType);
            end

            %% 4. Post-processing, Plotting and Printing
            h = obj.zLayer;
            plotMatlab = 1 && ~batchStartupOptionUsed;
            printParaview = true;
            kappa = 1;

            % Mid plane heights
            for i = 1:numel(h)-1
                z0 = h{i};
                z1 = h{i+1};
                zMidPlane{i} = 0.5*(z0+z1);
            end

            % Loop through solutions
            for iSol = 1:solutionsToProcess

                % ========== EXTRACT SOLUTION FOR CURRENT TIME/MODE ==========
                if strcmpi(problemType, 'FREE_VIBRATIONS')
                    fprintf('\n--- Processing Mode %d (%.4f Hz) ---\n', iSol, freq_Hz(iSol));
                    current_x = modes(:, iSol);
                    suffix = sprintf('_Mode_%d', iSol);

                elseif strcmpi(problemType, 'FORCED_VIBRATIONS')
                    time_idx = timeInstants(iSol);
                    current_time = time(time_idx);
                    fprintf('\n--- Processing Time Step %d / %d (t = %.4f s) ---\n', ...
                        time_idx, nt, current_time);
                    current_x = modes(:, time_idx);
                    suffix = sprintf('_Time_%d', time_idx);

                else  % STATIC
                    fprintf('\n--- Processing Static Solution ---\n');
                    current_x = modes(:, 1);
                    suffix = '_Static';
                end

                % Extract components
                uF = current_x(1:nU, 1);
                tF = current_x(nU+1 : nU+nTheta, 1);
                wF = current_x(nU+nTheta+1 : nU+nTheta+nW, 1);

                % Update Displacements (u)
                uT = zeros(obj.uFun.nDofs, 1);
                uT(dofFU, 1) = uF;
                uT = reshape(uT, obj.uFun.ndimf, [])';
                obj.uFun.setFValues(uT);

                % Update Transverse Displacements (w)
                wT = zeros(obj.wFun.nDofs, 1);
                wT(dofFW, 1) = wF;
                wT = reshape(wT, [], obj.wFun.ndimf);
                obj.wFun.setFValues(wT);

                % Update Rotations (theta)
                thetaT = zeros(obj.thetaFun.nDofs, 1);
                thetaT(dofFT, 1) = tF;
                thetaT = reshape(thetaT, obj.thetaFun.ndimf, [])';
                obj.thetaFun.setFValues(thetaT);

                % Max displacements
                [maxW, maxU, maxTheta, locationW, locationU, locationTheta] = obj.findMaxDisplacements();

                %% Strain and Stress Calculation
                [epsilons] = obj.createEpsilons();

                % Determine stress state
                stressState = 'PLANE_STRESS';
                [strainFun, stressFun] = obj.createStrainStressFunctions(obj.zLayer, epsilons, stressState);
                % obj.zLayer // zMidPlane

                vonMises = obj.computeVonMises(stressFun,stressState);

                % % PLOT STRESS DISTRIBUTION THROUGH THICKNESS
                % % Obtain max VonMises value on node and plot stresses
                % % through thickness
                % [maxVonMises, idxVM] = max(vonMises{end}.fValues);
                % locationVM = obj.mesh.coord(idxVM,:);
                % fprintf('Found maximum von Mises value on node %d (%.6e)\n', idxVM, maxVonMises);
                % fprintf('Location: (%.4f, %.4f)\n \n', locationVM(1), locationVM(2));
                % 
                % obj.plotStressDistributionThroughThickness(idxVM, epsilons, stressState);
                

                % ================================================
                %                 PLOT AND PRINT
                % ================================================
                if plotMatlab
                    if strcmpi(problemType, 'FREE_VIBRATIONS')
                        titleSuffix = sprintf(' - Mode %d (%.2f Hz)', iSol, freq_Hz(iSol));
                    elseif strcmpi(problemType, 'FORCED_VIBRATIONS')
                        titleSuffix = sprintf(' - t = %.4f s', current_time);
                    else
                        titleSuffix = '';
                    end

                    % ====== Displacements ======
                    % ===========================
                    obj.customPlot(obj.uFun, {['u_{x}' titleSuffix], ['u_{y}' titleSuffix]});
                    obj.customPlot(obj.wFun, {['w' titleSuffix]});
                    obj.customPlot(obj.thetaFun, {['\theta_{x}' titleSuffix], ['\theta_{y}' titleSuffix]});

                    % ======== Strains and Stresses ========
                    % ======================================
                    if strcmpi(stressState, 'PLANE_STRESS')
                        obj.customPlot(strainFun{kappa}, {'\epsilon_{xx}', '\epsilon_{yy}', '\epsilon_{yz}', '\epsilon_{xz}', '\epsilon_{xy}'});
                        obj.customPlot(stressFun{kappa}, {'\sigma_{xx}', '\sigma_{yy}', '\sigma_{yz}', '\sigma_{xz}', '\sigma_{xy}'});
                    else
                        obj.customPlot(strainFun{kappa}, {'\epsilon_{xx}', '\epsilon_{yy}', '\epsilon_{zz}', ...
                            '\epsilon_{yz}', '\epsilon_{xz}', '\epsilon_{xy}'});
                        obj.customPlot(stressFun{kappa}, {'\sigma_{xx}', '\sigma_{yy}', '\sigma_{zz}', ...
                            '\sigma_{yz}', '\sigma_{xz}', '\sigma_{xy}'});
                    end
                    % ======== Von Mises ========
                    % ===========================
                    obj.customPlot(vonMises{kappa},{'\sigma_{VM}'});
                end


                if printParaview
                    if strcmpi(problemType, 'STATIC')
                        % Static case: save in STATIC folder
                        baseFolder = 'STATIC';
                        if ~exist(baseFolder, 'dir')
                            mkdir(baseFolder);
                        end
                        outputPath = baseFolder;

                    else  % DYNAMIC CASES
                        % Dynamic case: DYNAMIC/problemType/TIMESTEP_X or MODE_X structure
                        baseDynamicFolder = 'DYNAMIC';
                        if ~exist(baseDynamicFolder, 'dir')
                            mkdir(baseDynamicFolder);
                        end

                        % Create subfolder for problem type
                        problemFolder = fullfile(baseDynamicFolder, problemType);
                        if ~exist(problemFolder, 'dir')
                            mkdir(problemFolder);
                        end

                        % Create subfolder for this solution
                        if strcmpi(problemType, 'FREE_VIBRATIONS')
                            solutionFolder = sprintf('MODE_%d', iSol);
                        elseif strcmpi(problemType, 'FORCED_VIBRATIONS')
                            solutionFolder = sprintf('TIMESTEP_%04d', time_idx);  % Zero-padded
                        else
                            solutionFolder = sprintf('SOLUTION_%d', iSol);
                        end

                        outputPath = fullfile(problemFolder, solutionFolder);
                        if ~exist(outputPath, 'dir')
                            mkdir(outputPath);
                        end
                    end

                    % ================================
                    % ========== SAVE FILES ==========
                    % ================================
                    fprintf('Saving Paraview files to: %s\n', outputPath);

                    obj.wFun.print(fullfile(outputPath, ['wfun_print' suffix]), 'Paraview');
                    obj.uFun.print(fullfile(outputPath, ['ufun_print' suffix]), 'Paraview');
                    obj.thetaFun.print(fullfile(outputPath, ['thetafun_print' suffix]), 'Paraview');
                    
                    stressFun{kappa}.print(fullfile(outputPath, ['stressfun' suffix]), 'Paraview')
                    strainFun{kappa}.print(fullfile(outputPath, ['strainfun' suffix]), 'Paraview')
                    vonMises{kappa}.print(fullfile(outputPath, ['VonMises' suffix]), 'Paraview')

                    % Save metadata
                    infoFile = fullfile(outputPath, 'info.txt');
                    fid = fopen(infoFile, 'w');
                    fprintf(fid, '===== ANALYSIS INFORMATION =====\n');
                    fprintf(fid, 'Problem Type: %s\n', problemType);
                    fprintf(fid, 'Date: %s\n', datestr(now));

                    if strcmpi(problemType, 'FORCED_VIBRATIONS')
                        fprintf(fid, '\n--- Time Step %d / %d ---\n', time_idx, nt);
                        fprintf(fid, 'Time: %.6f s\n', current_time);
                        fprintf(fid, 'Forcing frequency: %.4f Hz\n', f_force);
                        fprintf(fid, 'Max |w|: %.6e\n', maxW);
                        fprintf(fid, 'Max |u|: %.6e\n', maxU);
                        fprintf(fid, 'Max |θ|: %.6e\n', maxTheta);
                    elseif strcmpi(problemType, 'FREE_VIBRATIONS')
                        fprintf(fid, '\n--- Mode %d ---\n', iSol);
                        fprintf(fid, 'Natural frequency: %.4f Hz\n', freq_Hz(iSol));
                    end
                    fclose(fid);
                end
            end

        end

    end

    methods (Access = private)
        %% createMesh

        function createMesh(obj,meshtype)

            switch meshtype
                case 'unitTriangle'
                    obj.mesh = UnitTriangleMesh(50,50);
                case 'wingShape'

                    elements = 90;

                    fullmesh = TriangleMesh(18,10,elements,elements);
                    ls = obj.computeWingLevelSet(fullmesh);
                    sUm.backgroundMesh = fullmesh;
                    sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
                    uMesh              = UnfittedMesh(sUm);
                    uMesh.compute(ls);
                    wingMesh = uMesh.createInnerMesh();
                    obj.mesh = wingMesh;
            end
        end

        function ls = computeWingLevelSet(obj, mesh)
            gPar.type          = 'WingShape';
            gPar.xCoorCenter   = -0.01;
            gPar.yCoorCenter   = -0.01;
            gPar.chordRoot     = 7.3;
            gPar.chordTip      = 1.25;
            gPar.semiSpan      = 18.0;
            gPar.sweepDeg      = 25.0;
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsWing           = phiFun.fValues;
            ls = lsWing;
        end

        %% createMaterial
        function createMaterial(obj)
            % =========================================================================
            % MATERIAL PROPERTIES INPUT FORMAT
            % =========================================================================
            %
            % ISOTROPIC materials:
            %   E  = [E1; E2; ...; En]          % Young modulus per layer
            %   nu = [nu1; nu2; ...; nun]       % Poisson ratio per layer
            %   h  = [h1; h2; ...; hn]          % Thickness of each layer
            %   G is computed automatically as: G = E ./ (2*(1+nu))
            %
            % ORTHOTROPIC materials:
            %   E  = [E1 E2 E3; ...]            % Young moduli [longitudinal, transverse, through-thickness]
            %   nu = [nu12 nu13 nu23; ...]      % Poisson ratios
            %   G  = [G12 G13 G23; ...]         % Shear moduli
            %   h  = [h1; h2; ...; hn]          % Thickness of each layer
            %   Rotation = [theta1; theta2; ...] % Ply orientation in degrees
            %
            % =========================================================================
            % MATERIAL DATABASE
            % =========================================================================
            %
            % Available isotropic materials:
            %   'Aluminum' - Aluminum alloy (E=10.6 msi, nu=0.33)
            %   'Copper'   - Pure copper (E=18.0 msi, nu=0.33)
            %   'Steel'    - Structural steel (E=30.0 msi, nu=0.29)
            %
            % Available orthotropic composite materials:
            %   'AS'   - Graphite-Epoxy AS/3501 (High-strength carbon fiber)
            %   'EpT'  - Graphite-Epoxy T300/934 (Standard-modulus carbon fiber)
            %   'Ep1'  - Glass-Epoxy type 1 (E-glass, high strength)
            %   'Ep2'  - Glass-Epoxy type 2 (E-glass, lower modulus)
            %   'BrEp' - Boron-Epoxy (High stiffness, aerospace grade)
            %
            % =========================================================================

            % -------------------------------------------------------------------------
            % DATABASE MATERIALS
            % -------------------------------------------------------------------------

            materialCase = 'Composite';
            % 'Composite' // 'Aluminium'

            switch materialCase
                case 'Composite'
                    materialName = {'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT';
                        'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'};

                    max_thickness = 0.5;
                    Rotation = [0; 0; 45; -45; 90; 90; -45; 45; 0; 0];  % degrees
                    Rotation = 25*ones(size(materialName)) + Rotation;
                    
                    % Auto-distribute thickness
                    nLayers = length(materialName);
                    h = max_thickness / nLayers * ones(nLayers, 1);
                    
                    obj.dampingRatio = 0.015;
                    
                case 'Aluminium'
                    materialName = {'Aluminum'};
                    max_thickness = 0.5;
                    obj.dampingRatio = 0.01;
                   
                    % Auto-distribute thickness
                    nLayers = length(materialName);
                    h = max_thickness / nLayers * ones(nLayers, 1);
            end
             
            % materialName = {'Ep1'; 'Ep1'; 'Ep1'; 'Ep1'};
            % Rotation = [0;90;90;0];

            

            % Get material properties from database
            [E, nu, G, rho, type] = obj.getMaterialProperties(materialName);
            obj.materialProperties = type;

            % -------------------------------------------------------------------------
            % CUSTOM ORTHOTROPIC MATERIAL
            % -------------------------------------------------------------------------
            
            % obj.materialProperties = 'ISOTROPIC';
            % E  = [10; 3; 10];
            % nu = [0.3; 0.25; 0.3];
            % h  = [0.25; 0.5; 0.25];
            % G  = E ./ (2 .* (1 + nu));  % Auto-compute
            % Rotation = zeros(length(h), 1);  % Not used for isotropic

         
            % =========================================================================
            % STORE PROPERTIES IN OBJECT
            % =========================================================================

            obj.young   = ConstantFunction.create(E, obj.mesh);
            obj.poisson = ConstantFunction.create(nu, obj.mesh);
            obj.shear   = ConstantFunction.create(G, obj.mesh);
            obj.area    = ConstantFunction.create(h, obj.mesh);
            obj.inertia = ConstantFunction.create(1, obj.mesh);  
            obj.density = ConstantFunction.create(rho,obj.mesh);

            if strcmp(obj.materialProperties, 'ORTHOTROPIC')
                obj.rotation = Rotation;
            else
                obj.rotation = zeros(length(h), 1);
            end

            % =========================================================================
            % CREATE MATERIAL TENSOR
            % =========================================================================

            s.type    = 'MULTILAYER';
            s.ndim    = 2;  % Shell theory

            s.E       = E;
            s.nu      = nu;
            s.G       = G;
            s.h       = h;
            s.materialType = obj.materialProperties;

            if strcmp(obj.materialProperties, 'ORTHOTROPIC')
                s.rotation = Rotation;
            end

            tensor = Material.create(s);
            obj.material = tensor;

            % =========================================================================
            % COMPUTE ABDH TENSORS AND MATRIX 
            % =========================================================================

            nLayers = length(h);
            z_int = zeros(nLayers + 1, 1);
            z_int(1) = -sum(h)/2;
            for k = 1:nLayers
                z_int(k+1) = z_int(k) + h(k);
            end

            obj.zLayer = num2cell(z_int);   

            A_tens = zeros(3,3,3,3);
            B_tens = zeros(3,3,3,3);
            D_tens = zeros(3,3,3,3);
            H_tens = zeros(3,3,3,3);
            
            for k = 1:nLayers
                C_k = obj.material.getConstitutiveTensorForLayer(k);
                z0 = z_int(k);
                z1 = z_int(k+1);
                
                % Apply plane stress reduction (sigma_33 = 0) for A, B, D
                C_ps = obj.material.planeStressReduction(C_k);
                Qmatrix(:,:,k) = obj.material.tensorToVoigt(C_ps);
                
                A_tens = A_tens + C_ps * (z1 - z0);
                B_tens = B_tens + 0.5 * C_ps * (z1^2 - z0^2);
                D_tens = D_tens + 1/3 * C_ps * (z1^3 - z0^3);
                
                % H uses full 3D tensor (transverse shear, no plane stress reduction)
                H_tens = H_tens + C_k * (z1 - z0);
            end
            
            % Extract 2D in-plane components (a,b ∈ {1,2}) from the plane-stress condensed tensor
            A_2D = A_tens(1:2, 1:2, 1:2, 1:2);
            B_2D = B_tens(1:2, 1:2, 1:2, 1:2);
            D_2D = D_tens(1:2, 1:2, 1:2, 1:2);
            
            obj.A_tensor = ConstantFunction.create(A_2D, obj.mesh);
            obj.B_tensor = ConstantFunction.create(B_2D, obj.mesh);
            obj.D_tensor = ConstantFunction.create(D_2D, obj.mesh);
            
            % H: transverse shear (xz=13, yz=23) — no factor needed now that voigtToTensor
            % no longer divides shear components by 2
            H_2x2 = zeros(2,2);
            H_2x2(1,1) = H_tens(1,3,1,3);
            H_2x2(2,2) = H_tens(2,3,2,3);
            H_2x2(1,2) = H_tens(1,3,2,3);
            H_2x2(2,1) = H_tens(2,3,1,3);

            obj.H_tensor = ConstantFunction.create(H_2x2, obj.mesh);

            % Calculate Shear correction factor 
            index = [1,2,6];
            A_matrix = obj.material.tensorToVoigt(A_tens);
            B_matrix = obj.material.tensorToVoigt(B_tens);
            D_matrix = obj.material.tensorToVoigt(D_tens);
            F_bar = A_matrix(4,4) - A_matrix(4,5)^2 / A_matrix(5,5);

            ABD = [A_matrix(index,index), B_matrix(index,index); B_matrix(index,index), D_matrix(index,index)];
            Compliance = inv(ABD); 
            beta_matrix = Compliance(1:3,4:6);
            delta_matrix = Compliance(4:6,4:6);
            beta1i = beta_matrix(1,:);
            delta_1i = delta_matrix(1,:);
            H = zeros(1,nLayers);
            J = zeros(1,nLayers);

            for k = 1:nLayers
                H(k) = Qmatrix(1,1,k)*beta1i(1) + Qmatrix(1,2,k)*beta1i(2) + Qmatrix(1,6,k)*beta1i(3);
                J(k) = Qmatrix(1,1,k)*delta_1i(1) + Qmatrix(1,2,k)*delta_1i(2) + Qmatrix(1,6,k)*delta_1i(3);
            end
             
            sum_term = 0;
            for k = 1:nLayers
                zk  = obj.zLayer{k};
                zk1 = obj.zLayer{k+1};

                Tk = 0;
                Uk = 0;
                if k > 1
                    for m = 1:(k-1)
                        dz  = obj.zLayer{m+1} - obj.zLayer{m};
                        dz2 = obj.zLayer{m+1}^2 - obj.zLayer{m}^2;   
                        Tk = Tk + H(m) * dz;
                        Uk = Uk + (J(m)/2) * dz2;
                    end
                end

                Hk = H(k);
                Jk = J(k);

                Pk = Tk^2 + Hk^2*zk^2 - 2*Tk*Hk*zk + Uk^2 + (Jk^2*zk^4)/4 ...
                    - Uk*Jk*zk^2 + 2*Tk*Uk - Tk*Jk*zk^2 ...
                    - 2*Hk*Uk*zk + Hk*Jk*zk^3;  

                Rk = 2*Tk*Hk - 2*Hk^2*zk + 2*Hk*Uk - Hk*Jk*zk^2;
                Vk = Hk^2 - (Jk^2*zk^2)/2 + Uk*Jk + Tk*Jk - zk*Hk*Jk;
                Wk = Hk * Jk;
                Xk = (Jk^2) / 4;

                poly_int = Pk*(zk1 - zk) + (Rk/2)*(zk1^2 - zk^2) + ...
                    (Vk/3)*(zk1^3 - zk^3) + (Wk/4)*(zk1^4 - zk^4) + ...
                    (Xk/5)*(zk1^5 - zk^5);

                layer_shear_stiffness = Qmatrix(4,4,k) - (Qmatrix(4,5,k)^2) / Qmatrix(5,5,k);
                sum_term = sum_term + (1 / layer_shear_stiffness) * poly_int;
            end
            kappa = (1 / F_bar) / sum_term;
            obj.shearCorrectionFactor = kappa;

            
        end

        %% NewmarkMethod
        function [next_modalDisp, next_modalVel, next_modalAcc] = NewmarkMethod(obj, dt)
            
            alpha = 1/2; 
            gamma = 1/2;

            a1 = (1-alpha)*dt;
            a2 = alpha*dt; 
            a3 = 2/(gamma*dt^2);
            a4 = a3*dt; 
            a5 = (1-gamma)/gamma;

            b0 = a3; 
            b1 = a4; 
            b2 = a5; 
            b3 = a2*a3; 
            b4 = (a2*a4-1);
            b5 = (a2*a5-a1);

            A1 = inv(b0*obj.massModal + b3*obj.dampingModal + obj.stiffnessModal);
            A2 = A1*obj.massModal;
            A3 = A1*obj.dampingModal;

            next_modalDisp = @ (ini_modalDisp, ini_modalVel, ini_modalAcc, nextFmodal_t) A1 * nextFmodal_t + A2 * (b0*ini_modalDisp ...
                + b1*ini_modalVel + b2*ini_modalAcc) + A3 * (b3*ini_modalDisp + b4*ini_modalVel + b5*ini_modalAcc);
            next_modalAcc = @ (nextmodalDisp, ini_modalDisp, ini_modalVel, ini_modalAcc) a3 * nextmodalDisp - a3 * ini_modalDisp ...
                - a4 * ini_modalVel - a5 * ini_modalAcc;
            next_modalVel = @ (ini_modalVel, ini_modalAcc, nextmodalAcc) ini_modalVel + a1 * ini_modalAcc + a2 * nextmodalAcc;
        end

        %% createGeneralBoundaryConditions
        function bc = createGeneralBoundaryConditions(obj,direct)
            TOL = 1e-12;
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            isLeft  = @(coor)  abs(coor(:,1)-xMin)< TOL;
            isRight = @(coor)  abs(coor(:,1)-xMax)< TOL;
            isTop   = @(coor)  abs(coor(:,2)-yMax)< TOL;
            isBotom = @(coor)  abs(coor(:,2)-yMin)< TOL;


            % Boundary conditions Case
            obj.bcCase = 1;     % CHANGE THIS VALUE TO SELECT CASE: 1, 2, or 3

            switch obj.bcCase
                case 1
                    % CASE 1: Todo empotrado
                    % Todas las esquinas empotradas, carga q (1e5)
                    sDir{1}.domain    = @(coor) isLeft(coor) | isRight(coor) | isBotom(coor) | isTop(coor);
                    sDir{1}.direction = direct;
                    sDir{1}.value     = 0;
                    sDir{1}.ndim      = length(direct);

                    dirichletFun = [];
                    for i = 1:numel(sDir)
                        dir = DirichletCondition(obj.mesh, sDir{i});
                        dirichletFun = [dirichletFun, dir];
                    end
                    s.dirichletFun = dirichletFun;
                    s.pointloadFun = [];

                case 2
                    % CASE 2: Empotrada izquierda, carga q
                    % Lado izquierdo empotrado, carga q (50000) en toda la placa
                    sDir{1}.domain    = @(coor) isLeft(coor);
                    sDir{1}.direction = direct;
                    sDir{1}.value     = 0;
                    sDir{1}.ndim      = length(direct);

                    dirichletFun = [];
                    for i = 1:numel(sDir)
                        dir = DirichletCondition(obj.mesh, sDir{i});
                        dirichletFun = [dirichletFun, dir];
                    end
                    s.dirichletFun = dirichletFun;
                    s.pointloadFun = [];

                case 3
                    % CASE 3: Empotrada izquierda, pointload derecha
                    % Lado izquierdo empotrado, pointload en el nodo central derecho. Carga q = 0.
                    sDir{1}.domain    = @(coor) isLeft(coor);
                    sDir{1}.direction = direct;
                    sDir{1}.value     = 0;
                    sDir{1}.ndim      = length(direct);

                    dirichletFun = [];
                    for i = 1:numel(sDir)
                        dir = DirichletCondition(obj.mesh, sDir{i});
                        dirichletFun = [dirichletFun, dir];
                    end
                    s.dirichletFun = dirichletFun;

                    if length(direct) == 1
                        rightNodes = find(isRight(obj.mesh.coord));
                        if isempty(rightNodes)
                            error('No nodes found on the right boundary.');
                        end
                        idx = ceil(numel(rightNodes)/2);
                        singleNode = rightNodes(idx);

                        sPL{1}.domain    = @(coor) (1:size(coor,1))'==singleNode;
                        sPL{1}.direction = 2;
                        sPL{1}.value     = 10000;

                        pointloadFun = [];
                        for i = 1:numel(sPL)
                            pl = TractionLoad2(obj.mesh, sPL{i}, 'DIRAC');
                            pointloadFun = [pointloadFun, pl];
                        end
                        s.pointloadFun = pointloadFun;
                    else
                        s.pointloadFun = [];
                    end
            end

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);

        end

        %% createRHS
        function RHS = createRHS(obj)

            p = ConstantFunction.create([0 0],obj.mesh);
            m = ConstantFunction.create([0 0],obj.mesh);
            switch obj.bcCase
                case 1
                    q = ConstantFunction.create(1e5,obj.mesh);
                case 2
                    q = ConstantFunction.create(50000,obj.mesh);
                case 3
                    q = ConstantFunction.create(0,obj.mesh);
            end

            fu = @(v) DP(p,v);
            RHSu = IntegrateRHS(fu,obj.uFun,obj.mesh,'Domain',2);
            RHSu = obj.reduceVector(RHSu,obj.bcU);

            ftheta   = @(v) DP(m,v);
            RHStheta = IntegrateRHS(ftheta,obj.thetaFun,obj.mesh,'Domain',2);
            RHStheta = obj.reduceVector(RHStheta,obj.bcT);

            fw = @(v) q.*v;
            RHSw = IntegrateRHS(fw,obj.wFun,obj.mesh,'Domain',2);
            RHSw = obj.reduceVector(RHSw,obj.bcW);

            if obj.bcCase == 3
                obj.computeForces();
                RHSq = obj.RHSq;
                RHSq = obj.reduceVector(RHSq,obj.bcW);
                RHSw = RHSw + RHSq;
            end            

            RHS = [RHSu;RHStheta;RHSw];
        end


        %% createLHS
        function LHS = createLHS(obj)
            A = obj.A_tensor;
            f = @(u,v) DDP(SymGrad(v),DDP(A,SymGrad(u)));

            Ku = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Ku = reduceMatrix(obj,Ku,obj.bcU,obj.bcU);

            D = obj.D_tensor;
            f = @(u,v) DDP(SymGrad(v),DDP(D,SymGrad(u)));

            Ktheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Ktheta = obj.reduceMatrix(Ktheta,obj.bcT,obj.bcT);

            B = obj.B_tensor;
            f = @ (u,v) DDP(SymGrad(v),DDP(B,SymGrad(u)));

            Zut = IntegrateLHS(f,obj.uFun,obj.thetaFun,obj.mesh,'Domain',2);
            Zut = obj.reduceMatrix(Zut,obj.bcU,obj.bcT);


            H = obj.H_tensor;
            f = @(u,v) DP(v,DP(H,u));

            Mtheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Mtheta = obj.reduceMatrix(Mtheta,obj.bcT,obj.bcT);

            f = @(u,v) DP(v,DP(H,Grad(u)), 1, 1);

            Nthetaw = IntegrateLHS(f,obj.thetaFun,obj.wFun,obj.mesh,'Domain',2);
            Nthetaw = obj.reduceMatrix(Nthetaw,obj.bcT,obj.bcW);

            f = @(u,v) DP(Grad(v),DP(H,Grad(u)), 2, 1);

            Kw = IntegrateLHS(f,obj.wFun,obj.wFun,obj.mesh,'Domain',2);
            Kw = obj.reduceMatrix(Kw,obj.bcW,obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            beta = obj.shearCorrectionFactor;

            Ztu = Zut';
            Zuw = zeros(nU,nW);
            LHS = [Ku Zut Zuw; Ztu (Ktheta+beta*Mtheta) beta*Nthetaw; Zuw' beta*Nthetaw' beta*Kw];
        end

        %% createMassLHS
        function MLHS = createMassLHS(obj)

            zInterface = cell2mat(obj.zLayer);
            nLayers = length(zInterface);
            Arho = 0; Brho = 0; Drho = 0; 

            for k = 1:nLayers-1
                z0 = zInterface(k);
                z1 = zInterface(k+1);
                rho = obj.density.constant(k);
                
                Arho = Arho + rho * (z1 - z0);
                Brho = Brho + 0.5 * rho * (z1^2 - z0^2);
                Drho = Drho + 1/3 * rho * (z1^3 - z0^3);
            end

            f = @ (u,v) Arho*DP(v,u);
            Mu = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Mu = reduceMatrix(obj,Mu,obj.bcU,obj.bcU);

            f = @ (u,v) Brho*DP(v,u);
            Mut = IntegrateLHS(f,obj.uFun,obj.thetaFun,obj.mesh,'Domain',2);
            Mut = reduceMatrix(obj,Mut,obj.bcU,obj.bcT);

            f = @ (u,v) Drho*DP(v,u);
            Mt = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Mt = reduceMatrix(obj,Mt,obj.bcT,obj.bcT);

            Mtu = Mut.';

            f = @ (u,v) Arho.*v.*u;
            Mw = IntegrateLHS(f,obj.wFun,obj.wFun,obj.mesh,'Domain',2);
            Mw = reduceMatrix(obj,Mw,obj.bcW,obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            Zuw = zeros(nU,nW);
            Ztw = zeros(nTheta,nW);
            Zwu = Zuw.';
            Zwt = Ztw.';

            MLHS = [Mu, Mut, Zuw;
                Mtu, Mt, Ztw;
                Zwu, Zwt, Mw];

        end

        %% createStrainStressFunctions
        function [strainFun, stressFun] = createStrainStressFunctions(obj, z, epsilons, stressState)
            eps_u_11 = epsilons(:, 1);
            eps_u_22 = epsilons(:, 2);
            eps_u_12 = epsilons(:, 3);
            eps_theta_11 = epsilons(:, 4);
            eps_theta_22 = epsilons(:, 5);
            eps_theta_12 = epsilons(:, 6);
            dw_dx = epsilons(:, 7);
            dw_dy = epsilons(:, 8);

            % Transverse shear strains (tensorial: ε = γ/2)
            epsilon_xz = 0.5 * (dw_dx + obj.thetaFun.fValues(:, 1));
            epsilon_yz = 0.5 * (dw_dy + obj.thetaFun.fValues(:, 2));

            nNodes = size(obj.mesh.coord, 1);
            zInterfaces = [obj.zLayer{:}];

            % Check for interfaces
            thereIsInterface = false;
            obj.interfaceIndex = [];

            if numel(zInterfaces) > 2
                internalInterfaces = zInterfaces(2:end-1);
                tol = 1e-10;

                for k = 1:numel(z)
                    if any(abs(internalInterfaces - z{k}) <= tol)
                        thereIsInterface = true;
                        obj.interfaceIndex(end+1) = k;
                        % fprintf('Found interface at z index: %d (z = %.6f)\n', k, z{k});
                    end
                end
            end

            % Initialize storage
            if ~thereIsInterface
                stressTemp = cell(1, numel(z));
                strainTemp = cell(1, numel(z));
                fprintf('No internal interfaces detected.\n');
            else
                stressTemp = cell(2, numel(z));
                strainTemp = cell(1, numel(z));
            end

            % ========== DETERMINE STRESS STATE ==========

            if strcmp(stressState, 'PLANE_STRESS')
                nComponents = 5;  % [ε11, ε22, ε23, ε13, ε12] - without ε33
                fprintf('Using PLANE STRESS formulation (5 components)\n');
            else
                nComponents = 6;  % [ε11, ε22, ε33, ε23, ε13, ε12]
                fprintf('Using FULL 3D formulation (6 components)\n');
            end

            % ========== LOOP THROUGH Z-LEVELS ==========
            for i = 1:numel(z)
                z_k = z{i};

                % Build strain vector in Voigt notation
                if strcmp(stressState, 'PLANE_STRESS')
                    % [ε11, ε22, ε23, ε13, ε12] - 5 components (no ε33)
                    strain_voigt = zeros(nNodes, 5);
                    strain_voigt(:, 1) = eps_u_11 + z_k * eps_theta_11;       % ε_xx
                    strain_voigt(:, 2) = eps_u_22 + z_k * eps_theta_22;       % ε_yy
                    strain_voigt(:, 3) = epsilon_yz;                           % ε_yz
                    strain_voigt(:, 4) = epsilon_xz;                           % ε_xz
                    strain_voigt(:, 5) = eps_u_12 + z_k * eps_theta_12;       % ε_xy
                else
                    % [ε11, ε22, ε33, ε23, ε13, ε12] - 6 components (full 3D)
                    strain_voigt = zeros(nNodes, 6);
                    strain_voigt(:, 1) = eps_u_11 + z_k * eps_theta_11;       % ε_xx
                    strain_voigt(:, 2) = eps_u_22 + z_k * eps_theta_22;       % ε_yy
                    strain_voigt(:, 3) = 0;                                    % ε_zz = 0
                    strain_voigt(:, 4) = epsilon_yz;                           % ε_yz
                    strain_voigt(:, 5) = epsilon_xz;                           % ε_xz
                    strain_voigt(:, 6) = eps_u_12 + z_k * eps_theta_12;       % ε_xy
                end

                strainTemp{i} = LagrangianFunction.create(obj.mesh, nComponents, 'P1');
                strainTemp{i}.setFValues(strain_voigt);

                % ========== CHECK IF INTERFACE ==========
                isInternalInterface = false;
                if numel(zInterfaces) > 2
                    internalInterfaces = zInterfaces(2:end-1);
                    tol = 1e-10;
                    if any(abs(internalInterfaces - z_k) <= tol)
                        isInternalInterface = true;
                    end
                end

                % ========== COMPUTE STRESSES ==========
                if ~isInternalInterface
                    % Single layer
                    kLayer = find(z_k >= zInterfaces(1:end-1) & z_k <= zInterfaces(2:end), 1);

                    % Get 6x6 Voigt constitutive matrix
                    C_matrix_full = obj.material.createConstitutiveMatrixForLayer(kLayer);

                    % Apply rotation if needed
                    if strcmp(obj.material.materialType, 'ORTHOTROPIC') && obj.material.rotation(kLayer) ~= 0
                        theta = deg2rad(obj.material.rotation(kLayer));
                        C_matrix_full = obj.material.rotateConstitutiveMatrix(C_matrix_full, theta);
                    end

                    % Reduce to plane stress if needed
                    if strcmp(stressState, 'PLANE_STRESS')
                        idx = [1, 2, 4, 5, 6];
                        C_matrix = C_matrix_full(idx, idx) - ...
                                  (C_matrix_full(idx, 3) * C_matrix_full(3, idx)) / C_matrix_full(3, 3);  % 5x5 matrix
                        strain_voigt(:,3) = strain_voigt(:,3) * 2;
                        strain_voigt(:,4) = strain_voigt(:,4) * 2;
                        strain_voigt(:,5) = strain_voigt(:,5) * 2;
                    else
                        C_matrix = C_matrix_full;  % 6x6 matrix
                        strain_voigt(:,4) = strain_voigt(:,4) * 2;
                        strain_voigt(:,5) = strain_voigt(:,5) * 2;
                        strain_voigt(:,6) = strain_voigt(:,6) * 2;
                    end

                    % Compute stress: 
                    stress_voigt = (C_matrix * strain_voigt')';  % nNodes x nComponents

                    stressTemp{1, i} = LagrangianFunction.create(obj.mesh, nComponents, 'P1');
                    stressTemp{1, i}.setFValues(stress_voigt);

                else
                    % Interface: two layers
                    kLayerBottom = find(z_k >= zInterfaces(1:end-1) & z_k <= zInterfaces(2:end), 1);
                    kLayerTop = kLayerBottom + 1;

                    fprintf('  Interface at z = %.6f: Layer %d (top) and Layer %d (bottom)\n', ...
                        z_k, kLayerBottom, kLayerTop);

                    % Bottom layer
                    C_bottom_full = obj.material.createConstitutiveMatrixForLayer(kLayerBottom);
                    if strcmp(obj.material.materialType, 'ORTHOTROPIC') && obj.material.rotation(kLayerBottom) ~= 0
                        theta = deg2rad(obj.material.rotation(kLayerBottom));
                        C_bottom_full = obj.material.rotateConstitutiveMatrix(C_bottom_full, theta);
                    end

                    if strcmp(stressState, 'PLANE_STRESS')
                        idx = [1, 2, 4, 5, 6];
                        C_bottom = C_bottom_full(idx, idx) - ...
                            (C_bottom_full(idx, 3) * C_bottom_full(3, idx)) / C_bottom_full(3, 3);
                        strain_voigt(:,3) = strain_voigt(:,3) * 2;
                        strain_voigt(:,4) = strain_voigt(:,4) * 2;
                        strain_voigt(:,5) = strain_voigt(:,5) * 2;
                    else
                        C_bottom = C_bottom_full;
                        strain_voigt(:,4) = strain_voigt(:,4) * 2;
                        strain_voigt(:,5) = strain_voigt(:,5) * 2;
                        strain_voigt(:,6) = strain_voigt(:,6) * 2;
                    end
                    stress_bottom = (C_bottom * strain_voigt')';

                    stressTemp{1, i} = LagrangianFunction.create(obj.mesh, nComponents, 'P1');
                    stressTemp{1, i}.setFValues(stress_bottom);

                    % Top layer
                    C_top_full = obj.material.createConstitutiveMatrixForLayer(kLayerTop);
                    if strcmp(obj.material.materialType, 'ORTHOTROPIC') && obj.material.rotation(kLayerTop) ~= 0
                        theta = deg2rad(obj.material.rotation(kLayerTop));
                        C_top_full = obj.material.rotateConstitutiveMatrix(C_top_full, theta);
                    end

                    if strcmp(stressState, 'PLANE_STRESS')
                        idx = [1, 2, 4, 5, 6];
                        C_top = C_top_full(idx, idx) - ...
                            (C_top_full(idx, 3) * C_top_full(3, idx)) / C_top_full(3, 3);
                    else
                        C_top = C_top_full;
                    end
                    stress_top = (C_top * strain_voigt')';

                    stressTemp{2, i} = LagrangianFunction.create(obj.mesh, nComponents, 'P1');
                    stressTemp{2, i}.setFValues(stress_top);
                end
            end

            % ========== REORDER AND FLATTEN OUTPUT ==========
            count = 0;
            stressFun = {};
            strainFun = {};

            fprintf('\n===== OUTPUT ORGANIZATION =====\n');

            for i = 1:numel(z)
                z_k = z{i};
                isInterface = ~isempty(obj.interfaceIndex) && any(i == obj.interfaceIndex);

                if ~isInterface
                    count = count + 1;
                    strainFun{count} = strainTemp{i}; 
                    stressFun{count} = stressTemp{1, i}; 

                    kLayer = find(z_k >= zInterfaces(1:end-1) & z_k <= zInterfaces(2:end), 1);
                    fprintf('Position {%d}: z = %.6f, Layer %d\n', count, z_k, kLayer);
                else
                    kLayerBottom = find(z_k >= zInterfaces(1:end-1) & z_k <= zInterfaces(2:end), 1);
                    kLayerTop = kLayerBottom + 1;

                    % Top of bottom layer
                    count = count + 1;
                    strainFun{count} = strainTemp{i}; 
                    stressFun{count} = stressTemp{1, i}; 
                    fprintf('Position {%d}: z = %.6f, TOP of Layer %d\n', count, z_k, kLayerBottom);

                    % Bottom of top layer
                    count = count + 1;
                    strainFun{count} = strainTemp{i}; 
                    stressFun{count} = stressTemp{2, i}; 
                    fprintf('Position {%d}: z = %.6f, BOTTOM of Layer %d\n', count, z_k, kLayerTop);
                end
                
                
            end

            fprintf('\nTotal positions: %d\n', count);
            fprintf('==========================================\n\n');
        end

        %% internalForces
        function [Nfun, Mfun, Qfun] = internalForces(obj,epsilonU_nodal,epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal)

            Nvalues = obj.A_Matrix*epsilonU_nodal + obj.B_Matrix*epsilonTheta_nodal;
            Mvalues = obj.B_Matrix*epsilonU_nodal + obj.D_Matrix*epsilonTheta_nodal;
            Qvalues = obj.H_Matrix*[dw_dx_nodal.' ; dw_dy_nodal.'];

            Nfun = LagrangianFunction.create(obj.mesh,length(Nvalues(:,1)),'P1');
            Mfun = LagrangianFunction.create(obj.mesh,length(Mvalues(:,1)),'P1');
            Qfun = LagrangianFunction.create(obj.mesh,length(Qvalues(:,1)),'P1');
            Nfun.setFValues(Nvalues.');
            Mfun.setFValues(Mvalues.');
            Qfun.setFValues(Qvalues.');

        end


        %% createEpsilons
        function [epsilons] = createEpsilons(obj)
            GradSymU = SymGrad(obj.uFun);           % [ε_xx, ε_yy, ε_xy]
            GradSymTheta = SymGrad(obj.thetaFun);
            GradW = Grad(obj.wFun);                 % [∂w/∂x, ∂w/∂y]

            GradSymU = GradSymU.project('P1');
            GradSymTheta = GradSymTheta.project('P1');
            GradW = GradW.project('P1');


            % ========== ASSEMBLE OUTPUT MATRIX ==========
            
            epsilons(:,:, 1) = GradSymU.fValues(:,1);      % eps_u_11
            epsilons(:,:, 2) = GradSymU.fValues(:,4);      % eps_u_22
            epsilons(:,:, 3) = GradSymU.fValues(:,3);      % eps_u_12
            epsilons(:,:, 4) = GradSymTheta.fValues(:,1);  % eps_theta_11
            epsilons(:,:, 5) = GradSymTheta.fValues(:,4)';  % eps_theta_22
            epsilons(:,:, 6) = GradSymTheta.fValues(:,3)';  % eps_theta_12
            epsilons(:,:, 7) = GradW.fValues(:,1);                % dw_dx
            epsilons(:,:, 8) = GradW.fValues(:,2);                % dw_dy


        end

        %% createSolutionField
        function createSolutionField(obj)
           obj.uFun     = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.thetaFun = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.wFun     = LagrangianFunction.create(obj.mesh,1,'P1');
        end
        
        
        %% getMaterialProperties
        function [E, nu, G, density, type] = getMaterialProperties(obj,materialName)
            % ──────────────────────────────────────────────────────────────────────
            % AVAILABLE MATERIALS
            % ──────────────────────────────────────────────────────────────────────
            %
            % ISOTROPIC:
            %   'Aluminum', 'Copper', 'Steel'
            %
            % ORTHOTROPIC:
            %   'AS'   - Graphite-Epoxy AS/3501
            %   'EpT'  - Graphite-Epoxy T300/934
            %   'Ep1'  - Glass-Epoxy type 1
            %   'Ep2'  - Glass-Epoxy type 2
            %   'BrEp' - Boron-Epoxy

            % Conversion: 1 msi = 6894.76 MPa
            msi_to_Pa = 6.89476e9;

            % ==================== MATERIAL DATABASE ====================
            % ISOTROPIC
            db.Aluminum.type = 'ISOTROPIC';
            db.Aluminum.E  = 10.6 * msi_to_Pa;
            db.Aluminum.nu = 0.33;
            db.Aluminum.G  = 3.38 * msi_to_Pa;
            db.Aluminum.density = 2700;  % kg/m^3

            db.Copper.type = 'ISOTROPIC';
            db.Copper.E  = 18.0 * msi_to_Pa;
            db.Copper.nu = 0.33;
            db.Copper.G  = 6.39 * msi_to_Pa;
            db.Copper.density = 8960;

            db.Steel.type = 'ISOTROPIC';
            db.Steel.E  = 30.0 * msi_to_Pa;
            db.Steel.nu = 0.29;
            db.Steel.G  = 11.24 * msi_to_Pa;
            db.Steel.density = 7850;

            db.Al7075.type = 'ISOTROPIC';
            db.Al7075.E  = 67545e6;
            db.Al7075.nu = 0.33;
            db.Al7075.G  = 25393e6;
            db.Al7075.density = 2751;  % kg/m^3

            % ORTHOTROPIC
            db.AS.type = 'ORTHOTROPIC';
            db.AS.E  = [20.0, 1.3, 1.3] * msi_to_Pa;
            db.AS.nu = [0.30, 0.30, 0.49];
            db.AS.G  = [1.03, 1.03, 0.90] * msi_to_Pa;
            db.AS.density = 1600;

            db.EpT.type = 'ORTHOTROPIC';
            db.EpT.E  = [19.0, 1.5, 1.5] * msi_to_Pa;
            db.EpT.nu = [0.22, 0.22, 0.49];
            db.EpT.G  = [1.00, 0.90, 0.90] * msi_to_Pa;
            db.EpT.density = 1600;

            db.Ep1.type = 'ORTHOTROPIC';
            db.Ep1.E  = [7.8, 2.6, 2.6] * msi_to_Pa;
            db.Ep1.nu = [0.25, 0.25, 0.34];
            db.Ep1.G  = [1.30, 1.30, 0.50] * msi_to_Pa;
            db.Ep1.density = 1900; 

            db.Ep2.type = 'ORTHOTROPIC';
            db.Ep2.E  = [5.6, 1.2, 1.3] * msi_to_Pa;
            db.Ep2.nu = [0.26, 0.26, 0.34];
            db.Ep2.G  = [0.60, 0.60, 0.50] * msi_to_Pa;
            db.Ep2.density = 2000;

            db.BrEp.type = 'ORTHOTROPIC';
            db.BrEp.E  = [30.0, 3.0, 3.0] * msi_to_Pa;
            db.BrEp.nu = [0.30, 0.25, 0.25];
            db.BrEp.G  = [1.00, 1.00, 0.60] * msi_to_Pa;
            db.BrEp.density = 2000;

            db.T300_914_C.type = 'ORTHOTROPIC';
            db.T300_914_C.E  = [138.0, 11.0, 11.0] * 1e9;
            db.T300_914_C.nu = [0.28, 0.28, 0.40];
            db.T300_914_C.G  = [5.5, 5.5, 3.928] * 1e9;
            db.T300_914_C.density = 1580;

            % ==================== PROCESS INPUT ====================
            if ischar(materialName)
                materialName = {materialName};
                nMaterials = 1;
            elseif iscell(materialName)
                nMaterials = numel(materialName);
            else
                error('materialName must be char or cell array');
            end

            % Clean field names
            fieldNames = cell(nMaterials, 1);
            for i = 1:nMaterials
                fn = materialName{i};
                fn = strrep(fn, '-', '_');
                fn = strrep(fn, ' ', '_');
                fieldNames{i} = fn;
            end

            % ==================== CHECK MATERIALS EXIST ====================
            for i = 1:nMaterials
                fnames = fieldnames(db);
                found = 0;
                for j = 1:numel(fnames)
                    if strcmp(fieldNames{i}, fnames{j})
                        found = 1;
                        break;
                    end
                end
                if ~found
                    error('Material "%s" not found. Available: %s', ...
                        materialName{i}, sprintf('%s ', fnames{:}));
                end
            end

            % ==================== DETERMINE TYPE ====================
            hasIso = 0;
            hasOrtho = 0;
            for i = 1:nMaterials
                if strcmp(db.(fieldNames{i}).type, 'ISOTROPIC')
                    hasIso = 1;
                else
                    hasOrtho = 1;
                end
            end

            if hasOrtho
                type = 'ORTHOTROPIC';
                if hasIso
                    fprintf('Note: Mixed laminate. Converting isotropic to orthotropic.\n');
                end
            else
                type = 'ISOTROPIC';
            end

            % ==================== EXTRACT PROPERTIES ====================
            if strcmp(type, 'ISOTROPIC')
                E  = zeros(nMaterials, 1);
                nu = zeros(nMaterials, 1);
                G  = zeros(nMaterials, 1);
                density = zeros(nMaterials, 1);

                for i = 1:nMaterials
                    E(i)  = db.(fieldNames{i}).E;
                    nu(i) = db.(fieldNames{i}).nu;
                    G(i)  = db.(fieldNames{i}).G;
                    density(i) = db.(fieldNames{i}).density;
                end
            else
                E  = zeros(nMaterials, 3);
                nu = zeros(nMaterials, 3);
                G  = zeros(nMaterials, 3);
                density = zeros(nMaterials, 1);

                for i = 1:nMaterials
                    mat = db.(fieldNames{i});
                    if strcmp(mat.type, 'ISOTROPIC')
                        E(i,:)  = [mat.E, mat.E, mat.E];
                        nu(i,:) = [mat.nu, mat.nu, mat.nu];
                        G(i,:)  = [mat.G, mat.G, mat.G];
                        density(i) = db.(fieldNames{i}).density;
                    else
                        E(i,:)  = mat.E;
                        nu(i,:) = mat.nu;
                        G(i,:)  = mat.G;
                        density(i) = db.(fieldNames{i}).density;
                    end
                end
            end

            % Simplify single material output
            if nMaterials == 1
                s = size(E, 2);
                if s > 1
                    E = reshape(E, 1, s);
                    nu = reshape(nu, 1, s);
                    G = reshape(G, 1, s);
                    density = reshape(density,1,s);
                end
            end
        end

        %% computeVonMises 
        function vonMises = computeVonMises(obj,stressFun,stressState)
            ncells = numel(stressFun);

            for cell = 1:ncells
                stress = stressFun{cell}.fValues;
                if strcmp(stressState, 'PLANE_STRESS')
                    s1 = stress(:,1);
                    s2 = stress(:,2);
                    s12 = stress(:,3);
                    s23 = stress(:,4);
                    s31 = stress(:,5);

                    vonMises_vals = sqrt( ((s1 - s2).^2 + (s2).^2 + (s1).^2)/2 ...
                        + 3*(s12.^2 + s23.^2 + s31.^2) );
                else
                    s1 = stress(:,1);
                    s2 = stress(:,2);
                    s3 = stress(:,3);
                    s12 = stress(:,4);
                    s23 = stress(:,5);
                    s31 = stress(:,6);

                    vonMises_vals = sqrt( ((s1 - s2).^2 + (s2 - s3).^2 + (s3 - s1).^2)/2 ...
                        + 3*(s12.^2 + s23.^2 + s31.^2) );
                end

                vonMises{cell} = LagrangianFunction.create(obj.mesh,1,'P1');
                vonMises{cell}.setFValues(vonMises_vals);
            end

        end
        %% computeForces
        function computeForces(obj)
            bc  = obj.bcW;
            t   = bc.tractionFun;
            rhs = zeros(obj.wFun.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(obj.wFun);
                    rhs  = rhs + rhsi;
                end
            end
            if strcmp(obj.solverType,'REDUCED') 
                bc      = obj.bcW;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    % Convertir DOFs locales de w a globales (solo libres)
                    freeDofsW = obj.computeFreeDofs(obj.bcW);
                    globalDofsW = obj.localToGlobalDofs(freeDofsW, 'w');

                    R = -obj.lhs(globalDofsW, dirich)*dirichV;
                else
                    R = zeros(sum(obj.wFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
            obj.RHSq = rhs;
        end

        %% findMaxDisplacements
        function [maxW, maxU, maxTheta, locationW, locationU, locationTheta] = findMaxDisplacements(obj)
            % findMaxDisplacements - Find maximum displacements and their locations
            %
            % Outputs:
            %   maxW - Maximum absolute transverse displacement
            %   maxU - Maximum absolute in-plane displacement magnitude
            %   maxTheta - Maximum absolute rotation magnitude
            %   locationW - [x, y] coordinates of max w
            %   locationU - [x, y] coordinates of max |u|
            %   locationTheta - [x, y] coordinates of max |theta|

            % Get nodal coordinates
            coords = obj.mesh.coord;  % nNodes x 2 (x, y)

            % === TRANSVERSE DISPLACEMENT w ===
            w_values = obj.wFun.fValues;  % nNodes x 1
            [maxW, idxW] = max(abs(w_values));
            locationW = coords(idxW, :);

            % === IN-PLANE DISPLACEMENT u ===
            u_values = obj.uFun.fValues;  % nNodes x 2 (ux, uy)
            u_magnitude = sqrt(u_values(:, 1).^2 + u_values(:, 2).^2);
            [maxU, idxU] = max(u_magnitude);
            locationU = coords(idxU, :);

            % === ROTATIONS theta ===
            theta_values = obj.thetaFun.fValues;  % nNodes x 2 (theta_x, theta_y)
            theta_magnitude = sqrt(theta_values(:, 1).^2 + theta_values(:, 2).^2);
            [maxTheta, idxTheta] = max(theta_magnitude);
            locationTheta = coords(idxTheta, :);

            fprintf('\n===== MAXIMUM DISPLACEMENTS =====\n');
            fprintf('Max |w|:     %.6e at (%.4f, %.4f)\n', maxW, locationW(1), locationW(2));
            fprintf('Max |u|:     %.6e at (%.4f, %.4f)\n', maxU, locationU(1), locationU(2));
            fprintf('Max |theta|: %.6e at (%.4f, %.4f)\n', maxTheta, locationTheta(1), locationTheta(2));
            fprintf('=================================\n\n');
        end

        %% createBoundaryConditions
        function createBoundaryConditions(obj)            
            obj.bcU = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcT = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcW = obj.createGeneralBoundaryConditions([1]);            
        end
      

        %% localToGlobalDofs
        function globalDofs = localToGlobalDofs(obj, localDofs, field)
            % Convierte DOFs locales a DOFs globales en la matriz reducida
            % field puede ser 'u', 'theta' o 'w'

            nU = length(obj.bcU.free_dofs);
            nTheta = length(obj.bcT.free_dofs); 

            switch field
                case 'u'
                    offset = 0;
                case 'theta'
                    offset = nU;
                case 'w'
                    offset = nU + nTheta;
                otherwise
                    error('Field must be u, theta or w');
            end

            % Los localDofs ya son índices dentro del conjunto de DOFs libres
            % Solo necesitamos sumarles el offset
            globalDofs = length(localDofs) + offset;
        end
        
        
        %% computeFreeDofs
        function fD = computeFreeDofs(obj,bC)
            dofs = 1:bC.dirichletFun.nDofs;
            fD   = setdiff(dofs,bC.dirichlet_dofs); 
        end

        %% reduceMatrix
        function LHSred = reduceMatrix(obj,LHS,bcV,bcU)
            fdofV = obj.computeFreeDofs(bcV);
            fdofU = obj.computeFreeDofs(bcU);
            LHSred = LHS(fdofV,fdofU);
        end

        %% reduceVector
        function RHSred = reduceVector(obj,RHS,bc)
            fdofV = obj.computeFreeDofs(bc);
            RHSred = RHS(fdofV,1);
        end

         %% plotFigure
        % function customPlot(obj, fun, titles)
        %     % Inputs:
        %     %   fun    - LagrangianFunction to plot
        %     %   titles - Cell array of strings with titles for each subplot
        %     %
        %     % Example:
        %     %   customPlot(obj.uFun, {'u_{x}', 'u_{y}'});
        % 
        %     plot(fun);
        %     % colorbar;
        %     ax = findall(gcf, 'type', 'axes');
        %     ax = flipud(ax);
        % 
        %     for j = 1:length(ax)
        %         if j <= length(titles)
        %             title(ax(j), titles{j});
        %         end
        %     end
        % end

        function customPlot(obj, fun, titles)
            % Inputs:
            %   fun    - LagrangianFunction a plotear (escalar o vectorial)
            %   titles - Cell array con los títulos de cada subplot
            %
            % Ejemplo:
            %   customPlot(obj.sigmaVMFun, {'σ_{VM} - t = 1.2527 s'});

            plot(fun);                  % ← Se genera el plot primero
            drawnow;                    % ← Forzamos que se creen todos los objetos gráficos

            ax = findall(gcf, 'type', 'axes');
            ax = flipud(ax);

            for j = 1:length(ax)
                % Título
                if j <= length(titles)
                    title(ax(j), titles{j}, 'FontSize', 12, 'FontWeight', 'bold');
                end

                % === RECOLECTAR TODOS LOS DATOS DE COLOR (más robusto) ===
                cdata_all = [];

                % 1. Patch (el más común en mallas FEM como tu σ_VM)
                patches = findobj(ax(j), 'Type', 'Patch');
                for p = patches'
                    if isprop(p, 'FaceVertexCData') && ~isempty(p.FaceVertexCData)
                        cdata_all = [cdata_all; p.FaceVertexCData(:)];
                    end
                end

                % 2. Surface (trisurf, surf, etc.)
                surfaces = findobj(ax(j), 'Type', 'Surface');
                for s = surfaces'
                    if isprop(s, 'CData') && ~isempty(s.CData)
                        cdata_all = [cdata_all; s.CData(:)];
                    end
                end

                % 3. Image (por si acaso)
                images = findobj(ax(j), 'Type', 'Image');
                for im = images'
                    if isprop(im, 'CData') && ~isempty(im.CData)
                        cdata_all = [cdata_all; im.CData(:)];
                    end
                end

                % === AJUSTAR RANGO DE COLOR AL MÍNIMO Y MÁXIMO REAL ===
                if ~isempty(cdata_all)
                    cmin = min(cdata_all(:));
                    cmax = max(cdata_all(:));

                    if cmax > cmin + 1e-8
                        clim(ax(j), [cmin cmax]);     % rango exacto min → max
                    else
                        clim(ax(j), [cmin-1 cmax+1]); % evita rango degenerado
                    end
                end

                % % Colorbar opcional pero muy útil
                % cb = colorbar(ax(j));
                % cb.Label.String = titles{j};
                % cb.FontSize = 11;
            end
        end

        %% plotStressDistributionThroughThickness
        function plotStressDistributionThroughThickness(obj, nodeIdx, epsilons, stressState)
            if nargin < 4
                stressState = 'PLANE_STRESS';
            end

            % ========== VALIDATION ==========
            if nodeIdx > size(obj.mesh.coord, 1) || nodeIdx < 1
                error('Invalid node index: %d', nodeIdx);
            end

            % ========== GET COORDINATES OF NODE ==========
            coords = obj.mesh.coord;
            nodeCoords = coords(nodeIdx, :);

            % ========== GET LAYER INFORMATION ==========
            zInterfaces = [obj.zLayer{:}];
            nLayers = length(zInterfaces) - 1;

            % ========== BUILD UNIQUE z_values ==========
            pointsPerLayer = 4;
            z_unique = [];
            for k = 1:nLayers
                z_min = zInterfaces(k);
                z_max = zInterfaces(k+1);
                
                if k == 1
                    z_layer = linspace(z_min, z_max, pointsPerLayer)';
                else
                    z_layer = linspace(z_min, z_max, pointsPerLayer)';
                    z_layer = z_layer(2:end); 
                end
                
                z_unique = [z_unique; z_layer];
            end

            % ========== COMPUTE STRESSES ==========
            [~, stressFun] = obj.createStrainStressFunctions(num2cell(z_unique), epsilons, stressState);
            nPoints = numel(stressFun);
            sigmax = zeros(nPoints, 1);
            sigmay = zeros(nPoints, 1);

            for i = 1:nPoints
                sigmax(i) = stressFun{i}.fValues(nodeIdx, 1);
                sigmay(i) = stressFun{i}.fValues(nodeIdx, 2);
            end

            % ========== BUILD z_plot  ==========
            internalInterfaces = zInterfaces(2:end-1);   % vacío cuando nLayers = 1

            z_plot = [];
            for i = 1:length(z_unique)
                z = z_unique(i);
                if ~isempty(internalInterfaces) && any(abs(z - internalInterfaces) < 1e-10)
                    z_plot = [z_plot; z; z];
                else
                    z_plot = [z_plot; z];
                end
            end

            % ========== CREATE PLOT ==========
            figure;

            % --- SUBPLOT 1: sigmaxx vs z ---
            subplot(1, 2, 1);
            plot(sigmax/1e6, z_plot, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 6);
            grid on;
            xlabel('σ_{xx} (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('z (m)', 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('σ_{xx} Through Thickness - Node %d\n(%.4f, %.4f)', ...
                nodeIdx, nodeCoords(1), nodeCoords(2)), ...
                'FontSize', 11, 'FontWeight', 'bold');

            hold on;
            xline(0, 'k--', 'LineWidth', 1.1);
            if ~isempty(internalInterfaces)
                yline(internalInterfaces, 'k:', 'LineWidth', 1.0, 'Alpha', 0.65);
            end
            hold off;
            
            % --- SUBPLOT 2: σyy vs z ---
            subplot(1, 2, 2);
            plot(sigmay/1e6, z_plot, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 6);
            grid on;
            xlabel('σ_{yy} (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('z (m)', 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('σ_{yy} Through Thickness - Node %d\n(%.4f, %.4f)', ...
                nodeIdx, nodeCoords(1), nodeCoords(2)), ...
                'FontSize', 11, 'FontWeight', 'bold');

            hold on;
            xline(0, 'k--', 'LineWidth', 1.1);
            if ~isempty(internalInterfaces)
                yline(internalInterfaces, 'k:', 'LineWidth', 1.0, 'Alpha', 0.65);
            end
            hold off;
            
        end


        %% animateDynamicResponse
        function animateDynamicResponse(obj, modes, time, nU, nTheta, nW)
            % animateDynamicResponse - Animated plot of w displacement vs time

            % ========== FIND CORNER NODE ==========
            coords = obj.mesh.coord;  % All nodes in mesh
            [max_x, ~] = max(coords(:, 1));
            [max_y, ~] = max(coords(:, 2));
            distances = sqrt((coords(:, 1) - max_x).^2 + (coords(:, 2) - max_y).^2);
            [~, corner_node_global] = min(distances);  % Global node index

          

            fprintf('\n===== CORNER NODE INFORMATION =====\n');
            fprintf('Global node index: %d\n', corner_node_global);
            fprintf('Coordinates: (%.4f, %.4f)\n', coords(corner_node_global, 1), coords(corner_node_global, 2));

            % ========== CONVERT TO LOCAL w DOF INDEX ==========
            % Check if this node has a free w DOF
            dofFW = obj.computeFreeDofs(obj.bcW);

            % Find position of corner_node_global in the free DOFs list
            corner_node_local = find(dofFW == corner_node_global, 1);

            if isempty(corner_node_local)
                error('Corner node %d has no free w DOF (constrained by BC)', corner_node_global);
            end

            fprintf('Local w DOF index: %d (out of %d free w DOFs)\n', corner_node_local, nW);
            fprintf('===================================\n\n');

            % ========== EXTRACT W DISPLACEMENT HISTORY ==========
            % w DOFs start after u and theta
            w_start_idx = nU + nTheta + 1;
            w_end_idx = nU + nTheta + nW;

            % Extract all w displacements
            w_all_nodes = modes(w_start_idx:w_end_idx, :);  % nW x nt

            fprintf('w_all_nodes size: %d x %d\n', size(w_all_nodes, 1), size(w_all_nodes, 2));

            % Get w history for the corner node (using LOCAL index)
            w_corner = w_all_nodes(corner_node_local, :);  % 1 x nt

            % ========== CREATE ANIMATED PLOT ==========
            figure;
            h_plot = plot(time(1), w_corner(1), 'b-', 'LineWidth', 2);
            hold on;
            h_point = plot(time(1), w_corner(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            grid on;
            xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Displacement w (m)', 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('Dynamic Response at Corner Node %d (%.2f, %.2f)', ...
                corner_node_global, coords(corner_node_global, 1), coords(corner_node_global, 2)), ...
                'FontSize', 14, 'FontWeight', 'bold');
            xlim([time(1), time(end)]);

            % Set y-limits with some margin
            y_margin = 0.1 * max(abs(w_corner));
            ylim([min(w_corner) - y_margin, max(w_corner) + y_margin]);

            % Add horizontal line at y=0
            plot([time(1), time(end)], [0, 0], 'k--', 'LineWidth', 0.5);

            % ========== ANIMATE ==========
            fprintf('Starting animation...\n');
            skip = max(1, round(length(time) / 200));  % Show ~200 frames max

            for i = 1:skip:length(time)
                set(h_plot, 'XData', time(1:i), 'YData', w_corner(1:i));
                set(h_point, 'XData', time(i), 'YData', w_corner(i));
                title(sprintf('t = %.4f s | dt = %.6e m (Node %d)', ...
                    time(end), time(2)-time(1), corner_node_global));
                drawnow;
                pause(0.01);
            end

            hold off;

            fprintf('✓ Animation complete!\n');

            % ========== PRINT STATISTICS ==========
            [max_w, idx_max] = max(w_corner);
            [min_w, idx_min] = min(w_corner);

            fprintf('\n===== DISPLACEMENT STATISTICS =====\n');
            fprintf('Maximum: %.6e m at t = %.4f s\n', max_w, time(idx_max));
            fprintf('Minimum: %.6e m at t = %.4f s\n', min_w, time(idx_min));
            fprintf('Mean:    %.6e m\n', mean(w_corner));
            fprintf('RMS:     %.6e m\n', rms(w_corner));
            fprintf('Range:   %.6e m\n', max_w - min_w);
            fprintf('===================================\n\n');
        end

       
    end
    
end