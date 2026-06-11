classdef Al7075_WF110_ < handle

    properties (Access = private)
        mesh
        young
        poisson
        density
        area    
        thickness
        shear
        inertia
        uFun
        thetaFun
        wFun
        bcU,bcT,bcW
        h
        cost
        constraint
        material, rotation, shearCorrectionFactor
        zLayer
        A_tensor, B_tensor, D_tensor, H_tensor
        materialProperties, materialLayers
        meshType
    end

    methods (Access = public)
        %% Al7075_WF110_
        function obj = Al7075_WF110_(optimCase, forceSymmetry)

            if nargin < 1
                optimCase = 'STATIC'; % 'STATIC' / 'DYNAMIC'
            end
            if nargin < 2
                forceSymmetry = true;
            end

            clc; close all;

            % =========================================================================
            % INPUTS
            % =========================================================================
            obj.materialLayers = {'Al7075'; 'WF110'; 'Al7075'};

            thicknessRelation = 0.8;  % WF110 must be 80% of total thickness

            nLayers = length(obj.materialLayers);
            obj.rotation = zeros(nLayers,1);

            isSymmetric = obj.checkLaminateSymmetry(obj.rotation);

            % Initialization
            h_max = 0.8;
            h0    = h_max * ones(nLayers,1) / nLayers;

            obj.createMesh('wingShape')   % 'triangleMesh' // 'wingShape'
            obj.createSolutionField()

            % =========================================================================
            % TIMING: start global timer
            % =========================================================================
            tTotal = tic;

            % ---- Initial evaluation ----
            fprintf('Computing initial solution...\n');
            tPhase = tic;
            maxval_own = obj.staticProblem(h0);
            m0         = obj.computeMass(h0);
            fprintf('  Initial evaluation done in %.2f s\n', toc(tPhase));
            fprintf('  Initial maxval (Al7075/WF110): %g\n', maxval_own);

            % ---- Compute composite laminate maxval (from TutorialShellsOptim) ----
            fprintf('Computing composite laminate reference displacement...\n');
            tPhase2 = tic;
            maxval = obj.computeCompositeMaxval();
            fprintf('  Composite maxval done in %.2f s\n', toc(tPhase2));
            fprintf('  Reference maxval (composite): %g\n\n', maxval);

            % =========================================================================
            % OPTIMIZATION SETUP
            % =========================================================================
            switch optimCase

                case 'STATIC'
                    if isSymmetric && forceSymmetry
                        nHalf  = ceil(nLayers/2);
                        h0_opt = h0(1:nHalf);

                        I = eye(nHalf);
                        if mod(nLayers,2) == 0
                            T = [I; flipud(I)];
                        else
                            T = [I; flipud(I(1:end-1,:))];
                        end

                        s.ub = h_max * ones(nHalf,1);
                        s.lb = 0.01  * ones(nHalf,1);

                        cost.cF = @(h_half) obj.computeMass(T * h_half) / m0;
                        cost.gF = @(h_half) T' * obj.GradComputeMass() / m0;

                        constraint.cF{1} = @(h_half) obj.staticProblem(T * h_half) / maxval - 0.75;
                        constraint.gF{1} = @(h_half) T' * obj.GradStaticProblem(T * h_half, nLayers) / maxval;

                        % Constraint 2: h_WF110 / h_total = thicknessRelation  (EQUALITY)
                        % h_total = 2*h1 + h2  (symmetric: [h1, h2, h1])
                        constraint.cF{2} = @(h_half) h_half(2) / (2*h_half(1) + h_half(2)) - thicknessRelation;
                        constraint.gF{2} = @(h_half) [ ...
                            -2*h_half(2) / (2*h_half(1) + h_half(2))^2; ...
                             2*h_half(1) / (2*h_half(1) + h_half(2))^2 ];

                    else
                        h0_opt = h0;
                        s.ub   = h_max * ones(nLayers,1);
                        s.lb   = 0.01  * ones(nLayers,1);

                        cost.cF = @(h) obj.computeMass(h) / m0;
                        cost.gF = @(h) obj.GradComputeMass() / m0;

                        constraint.cF{1} = @(h) obj.staticProblem(h) / maxval - 0.75;
                        constraint.gF{1} = @(h) obj.GradStaticProblem(h, nLayers) / maxval;

                        % Constraint 2: h_WF110 / h_total = thicknessRelation  (EQUALITY)
                        % h_total = h1 + h2 + h3
                        constraint.cF{2} = @(h) h(2) / sum(h) - thicknessRelation;
                        constraint.gF{2} = @(h) [ ...
                            -h(2) / sum(h)^2; ...
                            (sum(h) - h(2)) / sum(h)^2; ...
                            -h(2) / sum(h)^2 ];
                    end

                case 'DYNAMIC'
                    fprintf('Computing initial dynamic solution (normalization)...\n');
                    tPhase = tic;
                    normalizef = obj.dynamicProblem(h0);
                    fprintf('  Done in %.2f s  (f1 = %.4f Hz)\n\n', toc(tPhase), normalizef);

                    if isSymmetric && forceSymmetry
                        nHalf  = ceil(nLayers/2);
                        h0_opt = h0(1:nHalf);

                        I = eye(nHalf);
                        if mod(nLayers,2) == 0
                            T = [I; flipud(I)];
                        else
                            T = [I; flipud(I(1:end-1,:))];
                        end

                        s.ub = h_max * ones(nHalf,1);
                        s.lb = 0.01  * ones(nHalf,1);

                        cost.cF = @(h_half) obj.computeMass(T * h_half) / m0;
                        cost.gF = @(h_half) T' * obj.GradComputeMass() / m0;

                        constraint.cF{1} = @(h_half) (-obj.dynamicProblem(T * h_half) + 10) / normalizef;
                        constraint.gF{1} = @(h_half) -T' * obj.GradDynamicProblem(T * h_half, nLayers) / normalizef;

                        % Constraint 2: h_WF110 / h_total = thicknessRelation  (EQUALITY)
                        constraint.cF{2} = @(h_half) h_half(2) / (2*h_half(1) + h_half(2)) - thicknessRelation;
                        constraint.gF{2} = @(h_half) [ ...
                            -2*h_half(2) / (2*h_half(1) + h_half(2))^2; ...
                             2*h_half(1) / (2*h_half(1) + h_half(2))^2 ];

                    else
                        h0_opt = h0;
                        s.ub   = h_max * ones(nLayers,1);
                        s.lb   = 0.01  * ones(nLayers,1);

                        cost.cF = @(h) obj.computeMass(h) / m0;
                        cost.gF = @(h) obj.GradComputeMass() / m0;

                        constraint.cF{1} = @(h) (-obj.dynamicProblem(h) + 10) / normalizef;
                        constraint.gF{1} = @(h) -obj.GradDynamicProblem(h, nLayers) / normalizef;

                        % Constraint 2: h_WF110 / h_total = thicknessRelation  (EQUALITY)
                        constraint.cF{2} = @(h) h(2) / sum(h) - thicknessRelation;
                        constraint.gF{2} = @(h) [ ...
                            -h(2) / sum(h)^2; ...
                            (sum(h) - h(2)) / sum(h)^2; ...
                            -h(2) / sum(h)^2 ];
                    end

                otherwise
                    error('Unknown optimCase: "%s". Use ''STATIC'' or ''DYNAMIC''.', optimCase);
            end

            % =========================================================================
            % SOLVER
            % =========================================================================
            s.type           = "fmincon";
            s.maxIter        = 20;
            s.constraintCase = {'INEQUALITY', 'EQUALITY'};
            s.tolerance      = 1e-5;

            cParams.cost         = cost;
            cParams.constraint   = constraint;
            cParams.initialGuess = h0_opt;
            cParams.settings     = s;
            cParams.printingPath = false;

            fprintf('Starting optimization (%s)...\n', optimCase);
            tOpt = tic;
            problem = AcademicProblem(cParams);
            problem.compute();
            tOptElapsed = toc(tOpt);
            fprintf('  Optimization finished in %.2f s\n\n', tOptElapsed);

            % =========================================================================
            % RECOVER FULL THICKNESS VECTOR
            % =========================================================================
            if isSymmetric && forceSymmetry
                h_vals = T * problem.result.fun.fValues;
            else
                h_vals = problem.result.fun.fValues;
            end

            % =========================================================================
            % TOTAL ELAPSED TIME
            % =========================================================================
            tTotalElapsed = toc(tTotal);
            [hh, mm, ss]  = obj.formatTime(tTotalElapsed);

            % =========================================================================
            % RESULTS REPORT
            % =========================================================================
            fprintf('\n==================================================\n');
            fprintf('               OPTIMIZATION RESULTS\n');
            fprintf('==================================================\n');

            for i = 1:length(h_vals)
                fprintf('  Layer %2d (%s): %10.6f m\n', i, obj.materialLayers{i}, h_vals(i));
            end

            fprintf('--------------------------------------------------\n');
            fprintf('  Total Thickness   : %10.6f m  (Limit: %.4f m)\n', sum(h_vals), h_max);
            fprintf('  WF110 fraction    : %10.4f %%  (Target: %.1f %%)\n', ...
                    h_vals(2)/sum(h_vals)*100, thicknessRelation*100);

            finalDisp = obj.staticProblem(h_vals);
            finalMass = obj.computeMass(h_vals);

            fprintf('  Final Disp (w)    : %10.6e m\n', finalDisp);
            fprintf('  Final Mass        : %10.6e kg  (Normalized: %.4f)\n', ...
                    finalMass, finalMass/m0);
            fprintf('--------------------------------------------------\n');
            fprintf('  Optim. time       : %10.2f s\n', tOptElapsed);
            fprintf('  TOTAL run time    : %02d h %02d min %05.2f s\n', hh, mm, ss);
            fprintf('==================================================\n\n');
        end
    end

    methods (Access = private)

        % =========================================================================
        % TIME FORMATTER
        % =========================================================================
        function [hh, mm, ss] = formatTime(~, t)
            hh = floor(t / 3600);
            mm = floor(mod(t, 3600) / 60);
            ss = mod(t, 60);
        end

        % =========================================================================
        % MESH
        % =========================================================================
        function createMesh(obj, meshtype)
            obj.meshType = meshtype;
            switch meshtype
                case 'triangleMesh'
                    obj.mesh = TriangleMesh(18,10,10,10);
                case 'wingShape'
                    elements = 50;
                    fullmesh = TriangleMesh(18,10,elements,elements);
                    ls = obj.computeWingLevelSet(fullmesh);
                    sUm.backgroundMesh = fullmesh;
                    sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
                    uMesh              = UnfittedMesh(sUm);
                    uMesh.compute(ls);
                    wingMesh = uMesh.createInnerMesh();
                    obj.mesh = wingMesh;
            end
            obj.area = obj.mesh.computeVolume();
        end

        function ls = computeWingLevelSet(obj, mesh)
            gPar.type        = 'WingShape';
            gPar.xCoorCenter = -0.01;
            gPar.yCoorCenter = -0.01;
            gPar.chordRoot   = 7.3;
            gPar.chordTip    = 1.25;
            gPar.semiSpan    = 18.0;
            gPar.sweepDeg    = 25.0;
            g                = GeometricalFunction(gPar);
            phiFun           = g.computeLevelSetFunction(mesh);
            ls               = phiFun.fValues;
        end

        % =========================================================================
        % SOLUTION FIELDS
        % =========================================================================
        function createSolutionField(obj)
            obj.uFun     = LagrangianFunction.create(obj.mesh, 2, 'P1');
            obj.thetaFun = LagrangianFunction.create(obj.mesh, 2, 'P1');
            obj.wFun     = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end

        % =========================================================================
        % STATIC PROBLEM
        % =========================================================================
        function wMax = staticProblem(obj, h)
            obj.createMaterialProperties(h);
            obj.createBoundaryConditions();
            LHS = obj.createLHS();
            RHS = obj.createRHS();

            x = LHS \ RHS;

            nU     = length(obj.computeFreeDofs(obj.bcU));
            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            wF    = x((nU+nTheta+1):(nU+nTheta+nW), 1);
            dofFW = obj.computeFreeDofs(obj.bcW);

            wT = zeros(obj.wFun.nDofs, 1);
            wT(dofFW, 1) = wF;
            wT = reshape(wT, [], obj.wFun.ndimf);
            obj.wFun.setFValues(wT);

            wMax = max(wT);
        end

        % =========================================================================
        % DYNAMIC PROBLEM
        % =========================================================================
        function freq1 = dynamicProblem(obj, h)
            obj.createMaterialProperties(h);
            obj.createBoundaryConditions();
            LHS  = obj.createLHS();
            MLHS = obj.createMassLHS();

            nModes = 3;
            [~, lambda]    = eigs(LHS, MLHS, nModes, 'smallestabs');
            omega_squared  = real(diag(lambda));
            omega_squared(omega_squared < 0) = 0;
            omega_n = sqrt(omega_squared);
            freq1   = omega_n(1) / (2*pi);
        end

        % =========================================================================
        % MASS
        % =========================================================================
        function mass = computeMass(obj, h)
            [~,~,~, rho, ~] = obj.getMaterialProperties(obj.materialLayers);
            mass = obj.area * sum(rho .* h);
        end

        function g = GradComputeMass(obj)
            [~,~,~, rho, ~] = obj.getMaterialProperties(obj.materialLayers);
            g = obj.area .* rho;
        end

        % =========================================================================
        % MATERIAL PROPERTIES
        % =========================================================================
        function createMaterialProperties(obj, h)
            materialName = obj.materialLayers;
            Rotation     = obj.rotation;

            [E, nu, G, rho, type] = obj.getMaterialProperties(materialName);
            obj.materialProperties = type;

            obj.young     = ConstantFunction.create(E,   obj.mesh);
            obj.poisson   = ConstantFunction.create(nu,  obj.mesh);
            obj.shear     = ConstantFunction.create(G,   obj.mesh);
            obj.thickness = ConstantFunction.create(h,   obj.mesh);
            obj.inertia   = ConstantFunction.create(1,   obj.mesh);
            obj.density   = ConstantFunction.create(rho, obj.mesh);

            s.type         = 'MULTILAYER';
            s.ndim         = 2;
            s.E            = E;
            s.nu           = nu;
            s.G            = G;
            s.h            = h;
            s.materialType = obj.materialProperties;
            if strcmp(obj.materialProperties, 'ORTHOTROPIC')
                s.rotation = Rotation;
            end
            tensor       = Material.create(s);
            obj.material = tensor;

            nLayers  = length(h);
            z_int    = zeros(nLayers+1, 1);
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
                C_k  = obj.material.getConstitutiveTensorForLayer(k);
                z0   = z_int(k);
                z1   = z_int(k+1);
                C_ps = obj.material.planeStressReduction(C_k);
                Qmatrix(:,:,k) = obj.material.tensorToVoigt(C_ps);
                A_tens = A_tens + C_ps * (z1 - z0);
                B_tens = B_tens + 0.5 * C_ps * (z1^2 - z0^2);
                D_tens = D_tens + (1/3) * C_ps * (z1^3 - z0^3);
                H_tens = H_tens + C_k  * (z1 - z0);
            end

            obj.A_tensor = ConstantFunction.create(A_tens(1:2,1:2,1:2,1:2), obj.mesh);
            obj.B_tensor = ConstantFunction.create(B_tens(1:2,1:2,1:2,1:2), obj.mesh);
            obj.D_tensor = ConstantFunction.create(D_tens(1:2,1:2,1:2,1:2), obj.mesh);

            H_2x2      = zeros(2,2);
            H_2x2(1,1) = H_tens(1,3,1,3);
            H_2x2(2,2) = H_tens(2,3,2,3);
            H_2x2(1,2) = H_tens(1,3,2,3);
            H_2x2(2,1) = H_tens(2,3,1,3);
            obj.H_tensor = ConstantFunction.create(H_2x2, obj.mesh);

            % Shear correction factor
            index    = [1,2,6];
            A_matrix = obj.material.tensorToVoigt(A_tens);
            B_matrix = obj.material.tensorToVoigt(B_tens);
            D_matrix = obj.material.tensorToVoigt(D_tens);
            F_bar    = A_matrix(4,4) - A_matrix(4,5)^2 / A_matrix(5,5);

            ABD        = [A_matrix(index,index), B_matrix(index,index); ...
                          B_matrix(index,index), D_matrix(index,index)];
            Compliance = inv(ABD);
            beta1i   = Compliance(1:3,4:6);    beta1i   = beta1i(1,:);
            delta_1i = Compliance(4:6,4:6);    delta_1i = delta_1i(1,:);

            H_scf = zeros(1,nLayers);
            J_scf = zeros(1,nLayers);
            for k = 1:nLayers
                H_scf(k) = Qmatrix(1,1,k)*beta1i(1) + Qmatrix(1,2,k)*beta1i(2) + Qmatrix(1,6,k)*beta1i(3);
                J_scf(k) = Qmatrix(1,1,k)*delta_1i(1) + Qmatrix(1,2,k)*delta_1i(2) + Qmatrix(1,6,k)*delta_1i(3);
            end

            sum_term = 0;
            for k = 1:nLayers
                zk  = obj.zLayer{k};
                zk1 = obj.zLayer{k+1};
                Tk  = 0; Uk = 0;
                if k > 1
                    for m = 1:(k-1)
                        dz  = obj.zLayer{m+1} - obj.zLayer{m};
                        dz2 = obj.zLayer{m+1}^2 - obj.zLayer{m}^2;
                        Tk  = Tk + H_scf(m) * dz;
                        Uk  = Uk + (J_scf(m)/2) * dz2;
                    end
                end
                Hk = H_scf(k); Jk = J_scf(k);
                Pk = Tk^2 + Hk^2*zk^2 - 2*Tk*Hk*zk + Uk^2 + (Jk^2*zk^4)/4 ...
                     - Uk*Jk*zk^2 + 2*Tk*Uk - Tk*Jk*zk^2 ...
                     - 2*Hk*Uk*zk + Hk*Jk*zk^3;
                Rk = 2*Tk*Hk - 2*Hk^2*zk + 2*Hk*Uk - Hk*Jk*zk^2;
                Vk = Hk^2 - (Jk^2*zk^2)/2 + Uk*Jk + Tk*Jk - zk*Hk*Jk;
                Wk = Hk * Jk;
                Xk = (Jk^2) / 4;
                poly_int = Pk*(zk1-zk) + (Rk/2)*(zk1^2-zk^2) + ...
                           (Vk/3)*(zk1^3-zk^3) + (Wk/4)*(zk1^4-zk^4) + ...
                           (Xk/5)*(zk1^5-zk^5);
                layer_shear_stiffness = Qmatrix(4,4,k) - Qmatrix(4,5,k)^2 / Qmatrix(5,5,k);
                sum_term = sum_term + (1/layer_shear_stiffness) * poly_int;
            end
            obj.shearCorrectionFactor = (1/F_bar) / sum_term;
        end

        % =========================================================================
        % LHS / RHS / BCs
        % =========================================================================
        function LHS = createLHS(obj)
            A = obj.A_tensor;
            f = @(u,v) DDP(SymGrad(v), DDP(A, SymGrad(u)));
            Ku = IntegrateLHS(f, obj.uFun, obj.uFun, obj.mesh, 'Domain', 2);
            Ku = obj.reduceMatrix(Ku, obj.bcU, obj.bcU);

            D = obj.D_tensor;
            f = @(u,v) DDP(SymGrad(v), DDP(D, SymGrad(u)));
            Ktheta = IntegrateLHS(f, obj.thetaFun, obj.thetaFun, obj.mesh, 'Domain', 2);
            Ktheta = obj.reduceMatrix(Ktheta, obj.bcT, obj.bcT);

            B = obj.B_tensor;
            f = @(u,v) DDP(SymGrad(v), DDP(B, SymGrad(u)));
            Zut = IntegrateLHS(f, obj.uFun, obj.thetaFun, obj.mesh, 'Domain', 2);
            Zut = obj.reduceMatrix(Zut, obj.bcU, obj.bcT);

            H = obj.H_tensor;
            f = @(u,v) DP(v, DP(H, u));
            Mtheta = IntegrateLHS(f, obj.thetaFun, obj.thetaFun, obj.mesh, 'Domain', 2);
            Mtheta = obj.reduceMatrix(Mtheta, obj.bcT, obj.bcT);

            f = @(u,v) DP(v, DP(H, Grad(u)), 1, 1);
            Nthetaw = IntegrateLHS(f, obj.thetaFun, obj.wFun, obj.mesh, 'Domain', 2);
            Nthetaw = obj.reduceMatrix(Nthetaw, obj.bcT, obj.bcW);

            f = @(u,v) DP(Grad(v), DP(H, Grad(u)), 2, 1);
            Kw = IntegrateLHS(f, obj.wFun, obj.wFun, obj.mesh, 'Domain', 2);
            Kw = obj.reduceMatrix(Kw, obj.bcW, obj.bcW);

            nU     = length(obj.computeFreeDofs(obj.bcU));
            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nW     = length(obj.computeFreeDofs(obj.bcW));
            beta   = obj.shearCorrectionFactor;

            Ztu = Zut';
            Zuw = zeros(nU, nW);
            LHS = [Ku,   Zut,                   Zuw; ...
                   Ztu, (Ktheta+beta*Mtheta),   beta*Nthetaw; ...
                   Zuw', beta*Nthetaw',          beta*Kw];
        end

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

            f = @(u,v) Arho * DP(v, u);
            Mu = IntegrateLHS(f, obj.uFun, obj.uFun, obj.mesh, 'Domain', 2);
            Mu = reduceMatrix(obj, Mu, obj.bcU, obj.bcU);

            f = @(u,v) Brho * DP(v, u);
            Mut = IntegrateLHS(f, obj.uFun, obj.thetaFun, obj.mesh, 'Domain', 2);
            Mut = reduceMatrix(obj, Mut, obj.bcU, obj.bcT);

            f = @(u,v) Drho * DP(v, u);
            Mt = IntegrateLHS(f, obj.thetaFun, obj.thetaFun, obj.mesh, 'Domain', 2);
            Mt = reduceMatrix(obj, Mt, obj.bcT, obj.bcT);

            Mtu = Mut.';

            f = @(u,v) Arho .* v .* u;
            Mw = IntegrateLHS(f, obj.wFun, obj.wFun, obj.mesh, 'Domain', 2);
            Mw = reduceMatrix(obj, Mw, obj.bcW, obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            Zuw = zeros(nU, nW);
            Ztw = zeros(nTheta, nW);
            Zwu = Zuw.';
            Zwt = Ztw.';

            MLHS = [Mu,   Mut,  Zuw; ...
                    Mtu,  Mt,   Ztw; ...
                    Zwu,  Zwt,  Mw];
        end

        function RHS = createRHS(obj)
            p = ConstantFunction.create([0 0], obj.mesh);
            m = ConstantFunction.create([0 0], obj.mesh);
            q = ConstantFunction.create(10240, obj.mesh);

            fu       = @(v) DP(p, v);
            RHSu     = IntegrateRHS(fu, obj.uFun, obj.mesh, 'Domain', 2);
            RHSu     = obj.reduceVector(RHSu, obj.bcU);

            ftheta   = @(v) DP(m, v);
            RHStheta = IntegrateRHS(ftheta, obj.thetaFun, obj.mesh, 'Domain', 2);
            RHStheta = obj.reduceVector(RHStheta, obj.bcT);

            fw       = @(v) q .* v;
            RHSw     = IntegrateRHS(fw, obj.wFun, obj.mesh, 'Domain', 2);
            RHSw     = obj.reduceVector(RHSw, obj.bcW);

            RHS = [RHSu; RHStheta; RHSw];
        end

        function createBoundaryConditions(obj)
            obj.bcU = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcT = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcW = obj.createGeneralBoundaryConditions([1]);
        end

        function bc = createGeneralBoundaryConditions(obj, direct)
            TOL  = 1e-12;
            xMin = min(obj.mesh.coord(:,1));
            isLeft = @(coor) abs(coor(:,1) - xMin) < TOL;

            sDir{1}.domain    = @(coor) isLeft(coor);
            sDir{1}.direction = direct;
            sDir{1}.value     = 0;
            sDir{1}.ndim      = length(direct);

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir          = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.periodicFun  = [];
            s.pointloadFun = [];
            s.mesh         = obj.mesh;
            bc             = BoundaryConditions(s);
        end

        function fD = computeFreeDofs(obj, bC)
            dofs = 1:bC.dirichletFun.nDofs;
            fD   = setdiff(dofs, bC.dirichlet_dofs);
        end

        function LHSred = reduceMatrix(obj, LHS, bcV, bcU)
            LHSred = LHS(obj.computeFreeDofs(bcV), obj.computeFreeDofs(bcU));
        end

        function RHSred = reduceVector(obj, RHS, bc)
            RHSred = RHS(obj.computeFreeDofs(bc), 1);
        end

        % =========================================================================
        % GRADIENTS (finite differences)
        % =========================================================================
        function g = GradStaticProblem(obj, h, nLayers)
            dh = 1e-6;
            g  = zeros(nLayers, 1);
            for i = 1:nLayers
                hf = h; hb = h;
                hf(i) = h(i) + dh;
                hb(i) = h(i) - dh;
                g(i) = (obj.staticProblem(hf) - obj.staticProblem(hb)) / (2*dh);
            end
        end

        function g = GradDynamicProblem(obj, h, nLayers)
            dh = 1e-6;
            g  = zeros(nLayers, 1);
            for i = 1:nLayers
                hf = h; hb = h;
                hf(i) = h(i) + dh;
                hb(i) = h(i) - dh;
                g(i) = (obj.dynamicProblem(hf) - obj.dynamicProblem(hb)) / (2*dh);
            end
        end

        % =========================================================================
        % MATERIAL DATABASE
        % =========================================================================
        function [E, nu, G, density, type] = getMaterialProperties(obj, materialName)
            msi_to_Pa = 6.89476e9;

            db.Aluminum.type    = 'ISOTROPIC';
            db.Aluminum.E       = 10.6 * msi_to_Pa;
            db.Aluminum.nu      = 0.33;
            db.Aluminum.G       = 3.38 * msi_to_Pa;
            db.Aluminum.density = 2700;

            db.Copper.type    = 'ISOTROPIC';
            db.Copper.E       = 18.0 * msi_to_Pa;
            db.Copper.nu      = 0.33;
            db.Copper.G       = 6.39 * msi_to_Pa;
            db.Copper.density = 8960;

            db.Steel.type    = 'ISOTROPIC';
            db.Steel.E       = 30.0 * msi_to_Pa;
            db.Steel.nu      = 0.29;
            db.Steel.G       = 11.24 * msi_to_Pa;
            db.Steel.density = 7850;

            db.Al7075.type    = 'ISOTROPIC';
            db.Al7075.E       = 67545e6;
            db.Al7075.nu      = 0.33;
            db.Al7075.G       = 25393e6;
            db.Al7075.density = 2751;

            db.WF110.type    = 'ISOTROPIC';
            db.WF110.E       = 194e6;
            db.WF110.G       = 67e6;
            db.WF110.nu      = db.WF110.E / (2*db.WF110.G) - 1;
            db.WF110.density = 110;

            db.AS.type    = 'ORTHOTROPIC';
            db.AS.E       = [20.0, 1.3, 1.3] * msi_to_Pa;
            db.AS.nu      = [0.30, 0.30, 0.49];
            db.AS.G       = [1.03, 1.03, 0.90] * msi_to_Pa;
            db.AS.density = 1600;

            db.EpT.type    = 'ORTHOTROPIC';
            db.EpT.E       = [19.0, 1.5, 1.5] * msi_to_Pa;
            db.EpT.nu      = [0.22, 0.22, 0.49];
            db.EpT.G       = [1.00, 0.90, 0.90] * msi_to_Pa;
            db.EpT.density = 1600;

            db.T300_914_C.type = 'ORTHOTROPIC';
            db.T300_914_C.E  = [138.0, 11.0, 11.0] * 1e9;
            db.T300_914_C.nu = [0.28, 0.28, 0.40];
            db.T300_914_C.G  = [5.5, 5.5, 3.928] * 1e9;
            db.T300_914_C.density = 1580;

            db.Ep1.type    = 'ORTHOTROPIC';
            db.Ep1.E       = [7.8, 2.6, 2.6] * msi_to_Pa;
            db.Ep1.nu      = [0.25, 0.25, 0.34];
            db.Ep1.G       = [1.30, 1.30, 0.50] * msi_to_Pa;
            db.Ep1.density = 1900;

            db.Ep2.type    = 'ORTHOTROPIC';
            db.Ep2.E       = [5.6, 1.2, 1.3] * msi_to_Pa;
            db.Ep2.nu      = [0.26, 0.26, 0.34];
            db.Ep2.G       = [0.60, 0.60, 0.50] * msi_to_Pa;
            db.Ep2.density = 2000;

            db.BrEp.type    = 'ORTHOTROPIC';
            db.BrEp.E       = [30.0, 3.0, 3.0] * msi_to_Pa;
            db.BrEp.nu      = [0.30, 0.25, 0.25];
            db.BrEp.G       = [1.00, 1.00, 0.60] * msi_to_Pa;
            db.BrEp.density = 2000;

            if ischar(materialName)
                materialName = {materialName};
            end
            nMaterials = numel(materialName);

            fieldNames = cell(nMaterials,1);
            for i = 1:nMaterials
                fn = strrep(materialName{i}, '-', '_');
                fn = strrep(fn, ' ', '_');
                fieldNames{i} = fn;
            end

            fnames = fieldnames(db);
            for i = 1:nMaterials
                if ~any(strcmp(fieldNames{i}, fnames))
                    error('Material "%s" not found. Available: %s', ...
                          materialName{i}, sprintf('%s ', fnames{:}));
                end
            end

            hasOrtho = any(cellfun(@(fn) strcmp(db.(fn).type,'ORTHOTROPIC'), fieldNames));
            hasIso   = any(cellfun(@(fn) strcmp(db.(fn).type,'ISOTROPIC'),   fieldNames));
            if hasOrtho
                type = 'ORTHOTROPIC';
                if hasIso
                    fprintf('Note: Mixed laminate. Converting isotropic to orthotropic.\n');
                end
            else
                type = 'ISOTROPIC';
            end

            if strcmp(type,'ISOTROPIC')
                E = zeros(nMaterials,1); nu = zeros(nMaterials,1);
                G = zeros(nMaterials,1); density = zeros(nMaterials,1);
                for i = 1:nMaterials
                    E(i)  = db.(fieldNames{i}).E;
                    nu(i) = db.(fieldNames{i}).nu;
                    G(i)  = db.(fieldNames{i}).G;
                    density(i) = db.(fieldNames{i}).density;
                end
            else
                E = zeros(nMaterials,3); nu = zeros(nMaterials,3);
                G = zeros(nMaterials,3); density = zeros(nMaterials,1);
                for i = 1:nMaterials
                    mat = db.(fieldNames{i});
                    if strcmp(mat.type,'ISOTROPIC')
                        E(i,:)  = [mat.E,  mat.E,  mat.E];
                        nu(i,:) = [mat.nu, mat.nu, mat.nu];
                        G(i,:)  = [mat.G,  mat.G,  mat.G];
                    else
                        E(i,:)  = mat.E;
                        nu(i,:) = mat.nu;
                        G(i,:)  = mat.G;
                    end
                    density(i) = mat.density;
                end
            end
        end

        % =========================================================================
        % SYMMETRY CHECK
        % =========================================================================
        function isSymmetric = checkLaminateSymmetry(obj, Rotation)
            nLayers     = length(obj.materialLayers);
            isSymmetric = true;
            tol         = 1e-6;
            for i = 1:floor(nLayers/2)
                j = nLayers - i + 1;
                if ~strcmp(obj.materialLayers{i}, obj.materialLayers{j})
                    isSymmetric = false; return;
                end
                if abs(Rotation(i) - Rotation(j)) > tol
                    isSymmetric = false; return;
                end
            end
        end

        % =========================================================================
        % COMPOSITE LAMINATE REFERENCE DISPLACEMENT
        % Replicates the TutorialShellsOptim configuration to compute the
        % initial displacement of the composite laminate (10 EpT layers).
        % =========================================================================
        function maxval = computeCompositeMaxval(obj)
            % Save current state
            savedMaterialLayers = obj.materialLayers;
            savedRotation       = obj.rotation;
            savedMaterialProps  = obj.materialProperties;

            % TutorialShellsOptim configuration
            compositeLayers = {'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'; ...
                               'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'};
            compositeRotation = [0; 0; 45; -45; 90; 90; -45; 45; 0; 0];
            compositeRotation = 25*ones(size(compositeLayers)) + compositeRotation;

            h_max_composite = 0.5;
            nLayersComposite = length(compositeLayers);
            h0_composite = h_max_composite * ones(nLayersComposite, 1) / nLayersComposite;

            % Temporarily set composite configuration
            obj.materialLayers = compositeLayers;
            obj.rotation       = compositeRotation;

            % Compute displacement
            maxval = obj.staticProblem(h0_composite);

            % Restore original state
            obj.materialLayers    = savedMaterialLayers;
            obj.rotation          = savedRotation;
            obj.materialProperties = savedMaterialProps;
        end

    end
end