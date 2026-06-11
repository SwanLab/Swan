classdef ValidationRefinement < handle


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
        lhs, Kw_full, RHSq
        solverType
        type, values, fun
        bcCase
        zLayer, interfaceIndex
        A_tensor, B_tensor, D_tensor, H_tensor
        materialProperties, stressCase
        material
        materialName
        naturalFrequencies
        modeShapes
        massModal, dampingModal, stiffnessModal, FModal
        dampingRatio
        shearCorrectionFactor
    end

    methods (Access = public)
        %% TutorialShells
        function obj = ValidationRefinement()
            close all; clc;

            mesh_sizes = [2, 10, 20, 30, 40, 50, 60, 70, 80];
            materials = {'Al7075', 'T300_914_C'};

            ref_w_Al = 4.26258377435646e-10;
            ref_w_Comp = 8.586112365969470e-10; % Calculated reference for composite

            errors_w = zeros(length(mesh_sizes), length(materials));
            values_w = zeros(length(mesh_sizes), length(materials));

            fig_err = figure('Color', 'w', 'Position', [200, 200, 800, 600]);
            hold on; 
            set(gca, 'YScale', 'log', 'FontSize', 16, 'LineWidth', 1.2, ...
                'TickLabelInterpreter', 'latex', 'XMinorTick', 'on', 'YMinorTick', 'on');
            grid on; grid minor;
            
            xlabel('Number of elements per side ($el$)', 'Interpreter', 'latex', 'FontSize', 18);
            ylabel('Relative Error in $w_{max}$', 'Interpreter', 'latex', 'FontSize', 18);
            title('\textbf{Mesh Refinement Convergence}', 'Interpreter', 'latex', 'FontSize', 22);
            
            colors = [0, 0.4470, 0.7410;      % Blue
                      0.8500, 0.3250, 0.0980]; % Orange
            markers = {'-o', '-s'};

            for m = 1:length(materials)
                obj.materialName = materials{m};
                if strcmp(obj.materialName, 'Al7075')
                    ref_w = ref_w_Al;
                else
                    ref_w = ref_w_Comp;
                end

                for idx_el = 1:length(mesh_sizes)
                    el = mesh_sizes(idx_el);

                    % 1. Initialization
                    obj.createMesh('unitTriangle', el)    % 'unitTriangle' // 'wingShape'
                    obj.createMaterial()
                    obj.createSolutionField()
                    obj.solverType = 'REDUCED';
                    obj.bcCase     = 1;


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

                    RHS = obj.createRHS();
                    x = LHS \ RHS;
                    current_x = x;

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

                    % Vertical displacement at the center
                    xc = 0.5; yc = 0.5;  % Para unitTriangle: a=1, b=1

                    coords = obj.mesh.coord;
                    dist   = sqrt((coords(:,1) - xc).^2 + (coords(:,2) - yc).^2);
                    [~, idx] = min(dist);

                    w_center = wT(idx);

                    % Strain and stresses calculation
                    [epsilons] = obj.createEpsilons();

                    % Determine stress state
                    stressState = 'PLANE_STRESS';
                    [strainFun, ~] = obj.createStrainStressFunctions(obj.zLayer, epsilons, stressState);

                    epsxx = strainFun{end}.fValues(idx,1);
                    epsyy = strainFun{end}.fValues(idx,2);

                    fprintf('Vertical displacement at center (%.3f, %.3f): \n  w = %.15e\n', ...
                        coords(idx,1), coords(idx,2), w_center);

                    % Imprimir deformaciones en el punto central (epsxx, epsyy, gxy, gxz, gyz)
                    % fprintf('Strains at center (coord [%.3f, %.3f]):\n', coords(idx,1), coords(idx,2));
                    fprintf('  eps_xx = %.15e at (%.3f, %.3f)\n', epsxx, coords(idx,1), coords(idx,2));
                    fprintf('  eps_yy = %.15e at (%.3f, %.3f)\n', epsyy, coords(idx,1), coords(idx,2));


                    coords = obj.mesh.coord;
                    [~, idx_00]  = min(sqrt((coords(:,1)-0  ).^2 + (coords(:,2)-0  ).^2));
                    [~, idx_0b2] = min(sqrt((coords(:,1)-0  ).^2 + (coords(:,2)-0.5).^2));
                    [~, idx_a20] = min(sqrt((coords(:,1)-0.5).^2 + (coords(:,2)-0  ).^2));

                    gxy_00  = 2 * strainFun{end}.fValues(idx_00,  5);
                    gxz_0b2 = 2 * strainFun{end}.fValues(idx_0b2, 4);
                    gyz_a20 = 2 * strainFun{end}.fValues(idx_a20, 3);

                    fprintf('  gxy  at (0, 0) [node (%.4f,%.4f)]: %+.15E\n', coords(idx_00, 1),  coords(idx_00, 2),  gxy_00);
                    fprintf('  gxz  at (0, b/2) [node (%.4f,%.4f)]: %+.15E\n', coords(idx_0b2,1),  coords(idx_0b2,2),  gxz_0b2);
                    fprintf('  gyz  at (a/2, 0) [node (%.4f,%.4f)]: %+.15E\n', coords(idx_a20,1),  coords(idx_a20,2),  gyz_a20);

                    % ─── Valores de referencia (analytical/reference) ───
                    if strcmp(obj.materialName, 'Al7075')
                        ref_w       =  4.26258377435646e-10;
                    elseif strcmp(obj.materialName, 'T300_914_C')
                        ref_w       =  8.586112365969470e-10;
                        
                    else
                        % Default/dummy values if another material is used
                        ref_w       = 0;
                    end

                    % ─── Errores relativos ───
                    err_w       = abs(w_center   - ref_w      ) / abs(ref_w      );

                                     
                    errors_w(idx_el, m) = err_w;
                    values_w(idx_el, m) = abs(w_center);
                end % End element loop

                figure(fig_err);
                plot(mesh_sizes, errors_w(:, m), markers{m}, 'LineWidth', 2, ...
                     'MarkerSize', 8, 'MarkerFaceColor', colors(m,:), ...
                     'Color', colors(m,:), ...
                     'DisplayName', sprintf('%s (Ref $w$: %.2e)', materials{m}, ref_w));
            end % End material loop

            figure(fig_err);
            lgd = legend('show', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best');
            lgd.Box = 'off';

            fig_w = figure('Color', 'w', 'Position', [250, 250, 800, 600]);
            hold on;
            set(gca, 'FontSize', 16, 'LineWidth', 1.2, ...
                'TickLabelInterpreter', 'latex', 'XMinorTick', 'on', 'YMinorTick', 'on');
            grid on; grid minor;
            xlabel('Number of elements per side ($el$)', 'Interpreter', 'latex', 'FontSize', 18);
            ylabel('Absolute value of $w_{max}$ [m]', 'Interpreter', 'latex', 'FontSize', 18);
            title('\textbf{Transverse Displacement Convergence}', 'Interpreter', 'latex', 'FontSize', 22);
            
            for m = 1:length(materials)
                if strcmp(materials{m}, 'Al7075')
                    ref_w = ref_w_Al;
                else
                    ref_w = ref_w_Comp;
                end
                
                plot(mesh_sizes, values_w(:, m), markers{m}, 'LineWidth', 2, ...
                     'MarkerSize', 8, 'MarkerFaceColor', colors(m,:), ...
                     'Color', colors(m,:), ...
                     'DisplayName', sprintf('%s (Calculated)', materials{m}));
                     
                yline(abs(ref_w), '--', 'Color', colors(m,:), 'LineWidth', 1.5, ...
                      'DisplayName', sprintf('%s (Analytical)', materials{m}));
            end
            
            lgd2 = legend('show', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best');
            lgd2.Box = 'off';
        end

    end

    methods (Access = private)
        %% createMesh

        function createMesh(obj,meshtype,el)

            if nargin < 3
                el = 80;
            end
            switch meshtype
                case 'unitTriangle'
                    obj.mesh = UnitTriangleMesh(el,el);
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
            %   'Al7075'
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

            materialName = {obj.materialName};
            max_thickness = 0.1;
            obj.dampingRatio = 0.01;

            nLayers = length(materialName);
            h = max_thickness / nLayers * ones(nLayers, 1);
                                   
            % Get material properties from database
            [E, nu, G, rho, type] = obj.getMaterialProperties(materialName);
            obj.materialProperties = type;
            Rotation = zeros(nLayers, 1);

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

            if length(direct) > 1
                % u_x / theta_x = 0 on top & bottom 
                sDir{1}.domain    = @(coor) isTop(coor) | isBotom(coor);
                sDir{1}.direction = direct(1);
                sDir{1}.value     = 0;
                sDir{1}.ndim      = length(direct);

                % u_y / theta_y = 0 on left & right 
                sDir{2}.domain    = @(coor) isLeft(coor) | isRight(coor);
                sDir{2}.direction = direct(2);
                sDir{2}.value     = 0;
                sDir{2}.ndim      = length(direct);
            else 
                % w = 0 on all boundaries 
                sDir{1}.domain    = @(coor) isTop(coor) | isBotom(coor) | isLeft(coor) | isRight(coor);
                sDir{1}.direction = direct;
                sDir{1}.value     = 0;
                sDir{1}.ndim      = length(direct);
            end

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            s.pointloadFun = [];
            if length(direct) == 1 && obj.bcCase == 2
                q0 = 1;
                Lx = 1; Ly = 1;
                
                % Apply sinusoidal load as nodal point loads
                nodal_q = q0 .* sin(pi .* obj.mesh.coord(:,1) ./ Lx) .* sin(pi .* obj.mesh.coord(:,2) ./ Ly);
                
                % Approximate integration (force = pressure * area)
                area_per_node = (Lx * Ly) / size(obj.mesh.coord,1);
                
                sPL{1}.domain    = @(coor) true(size(coor,1), 1);
                sPL{1}.direction = 1;
                sPL{1}.value     = nodal_q .* area_per_node;
                
                pointloadFun = [];
                for i = 1:numel(sPL)
                    pl = TractionLoad2(obj.mesh, sPL{i}, 'DIRAC');
                    pointloadFun = [pointloadFun, pl];
                end
                s.pointloadFun = pointloadFun;
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
                    q0 = 1;
                    Lx = 1;
                    Ly = 1;
                    fHandle = @(x) q0 .* sin(pi .* x(1,:,:) ./ Lx) .* sin(pi .* x(2,:,:) ./ Ly);
                    q = AnalyticalFunction.create(fHandle, obj.mesh);
                case 2
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

            obj.computeForces();
            RHSq = obj.RHSq;
            RHSq = obj.reduceVector(RHSq,obj.bcW);
            RHSw = RHSw + RHSq;

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
            obj.Kw_full = Kw;
            Kw = obj.reduceMatrix(Kw,obj.bcW,obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            beta = obj.shearCorrectionFactor;

            Ztu = Zut';
            Zuw = zeros(nU,nW);
            LHS = [Ku Zut Zuw; Ztu (Ktheta+beta*Mtheta) beta*Nthetaw; Zuw' beta*Nthetaw' beta*Kw];
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
                    % density is a scalar, no need to reshape
                end
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
                    freeDofsW = obj.computeFreeDofs(obj.bcW);
                    beta = obj.shearCorrectionFactor;
                    R = -beta * obj.Kw_full(freeDofsW, dirich) * dirichV;
                    rhs(freeDofsW) = rhs(freeDofsW) + R;
                end
            end
            obj.RHSq = rhs;
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
            globalDofs = localDofs + offset;
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

       
       
       
    end
end
    
