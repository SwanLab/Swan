classdef DifferentThetasTutorialShellsOptim_ < handle

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
        cost
        constraint
        material, rotation, shearCorrectionFactor
        zLayer
        A_tensor, B_tensor, D_tensor, H_tensor
        materialProperties, materialLayers
        meshType

        last_h_eval = []
        last_f1_eval = []
    end

    methods (Access = public)
        function obj = DifferentThetasTutorialShellsOptim_(optimCase_in, nOuter_in, deltaTheta0_in)

            clc; close all;
            % =========================================================================
            % TIMING: start global timer
            % =========================================================================
            tTotal = tic;

            if nargin > 0
                optimCase = optimCase_in;
            else
                optimCase = 'STATIC';
            end
            if nargin > 1
                nOuter = nOuter_in;
            else
                nOuter = 5;
            end
            if nargin > 2
                deltaTheta0 = deltaTheta0_in;
            else
                deltaTheta0 = 70;
            end

            obj.materialLayers = {'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C';
                                  'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'};

            Rotation = [0; 0; 45; -45; 90; 90; -45; 45; 0; 0];
            Rotation = 25*ones(size(obj.materialLayers)) + Rotation;
            obj.rotation = Rotation;

            forceSymmetry = true;
            isSymmetric   = obj.checkLaminateSymmetry(obj.rotation);

            h_max   = 0.5;
            nLayers = length(obj.materialLayers);
            h0      = h_max * ones(nLayers,1) / nLayers;
            % deltaTheta0 set from input argument

            obj.createMesh('wingShape')
            obj.createSolutionField()

            minDeltaTheta = -89.9;
            maxDeltaTheta =  89.9;

            normalizef = 20;

            % ---- Initial evaluation ----
            fprintf('Computing initial solution...\n');
            tPhase = tic;
            maxval = obj.staticProblem(h0, zeros(nLayers,1));
            m0     = obj.computeMass(h0);
            fprintf('  Initial evaluation done in %.2f s\n', toc(tPhase));
            fprintf('  Initial maxval: %g\n\n', maxval);

            switch optimCase
                case 'STATIC'
                    if isSymmetric && forceSymmetry
                        nHalf = ceil(nLayers/2);
                        I = eye(nHalf);
                        if mod(nLayers,2) == 0
                            T = [I; flipud(I)];
                        else
                            T = [I; flipud(I(1:end-1,:))];
                        end

                        % S: (nHalf x nOuter)
                        S = [eye(nOuter); zeros(nHalf-nOuter, nOuter)];

                        % x = [h(1:nHalf); deltaTheta(1:nOuter)]
                        x0_opt = [h0(1:nHalf); deltaTheta0*ones(nOuter,1)];

                        s.ub = [h_max        * ones(nHalf,1);  maxDeltaTheta * ones(nOuter,1)];
                        s.lb = [0.01         * ones(nHalf,1);  minDeltaTheta * ones(nOuter,1)];

                        cost.cF = @(x) obj.computeMass(T * x(1:nHalf)) / m0;
                        cost.gF = @(x) [T' * obj.GradComputeMass() / m0; zeros(nOuter,1)];

                        constraint.cF{1} = @(x) ...
                            obj.staticProblem(T * x(1:nHalf), T * S * x(nHalf+1:nHalf+nOuter)) / maxval - 0.75;
                        constraint.gF{1} = @(x) ...
                            obj.GradStaticProblem(x, nHalf, nOuter, T, S) / maxval;

                    else
                        S = [eye(nOuter); zeros(nLayers-nOuter, nOuter)];

                        x0_opt = [h0; deltaTheta0*ones(nOuter,1)];

                        s.ub = [h_max * ones(nLayers,1);  maxDeltaTheta * ones(nOuter,1)];
                        s.lb = [0.01  * ones(nLayers,1);  minDeltaTheta * ones(nOuter,1)];

                        cost.cF = @(x) obj.computeMass(x(1:nLayers)) / m0;
                        cost.gF = @(x) [obj.GradComputeMass() / m0; zeros(nOuter,1)];

                        constraint.cF{1} = @(x) ...
                            obj.staticProblem(x(1:nLayers), S * x(nLayers+1:nLayers+nOuter)) / maxval - 0.75;
                        constraint.gF{1} = @(x) ...
                            obj.GradStaticProblem(x, nLayers, nOuter, eye(nLayers), S) / maxval;
                    end

                case 'DYNAMIC'
                    if isSymmetric && forceSymmetry
                        nHalf = ceil(nLayers/2);
                        I = eye(nHalf);
                        if mod(nLayers,2) == 0
                            T = [I; flipud(I)];
                        else
                            T = [I; flipud(I(1:end-1,:))];
                        end

                        S = [eye(nOuter); zeros(nHalf-nOuter, nOuter)];

                        x0_opt = [h0(1:nHalf); deltaTheta0*ones(nOuter,1)];

                        s.ub = [h_max * ones(nHalf,1);  maxDeltaTheta * ones(nOuter,1)];
                        s.lb = [0.01  * ones(nHalf,1);  minDeltaTheta * ones(nOuter,1)];

                        cost.cF = @(x) obj.computeMass(T * x(1:nHalf)) / m0;
                        cost.gF = @(x) [T' * obj.GradComputeMass() / m0; zeros(nOuter,1)];

                        constraint.cF{1} = @(x) ...
                           - (obj.dynamicProblem(T * x(1:nHalf), T * S * x(nHalf+1:nHalf+nOuter)) - 10) / normalizef;
                        constraint.gF{1} = @(x) ...
                           - obj.GradDynamicProblem(x, nHalf, nOuter, T, S) / normalizef;

                        % constraint.cF{2} = @(x) -constraint.cF{1}(x) - 10/normalizef;
                        % constraint.gF{2} = @(x) -constraint.gF{1}(x);

                    else
                        S = [eye(nOuter); zeros(nLayers-nOuter, nOuter)];

                        x0_opt = [h0; deltaTheta0*ones(nOuter,1)];

                        s.ub = [h_max * ones(nLayers,1);  maxDeltaTheta * ones(nOuter,1)];
                        s.lb = [0.01  * ones(nLayers,1);  minDeltaTheta * ones(nOuter,1)];

                        cost.cF = @(x) obj.computeMass(x(1:nLayers)) / m0;
                        cost.gF = @(x) [obj.GradComputeMass() / m0; zeros(nOuter,1)];

                        constraint.cF{1} = @(x) ...
                            (obj.dynamicProblem(x(1:nLayers), S * x(nLayers+1:nLayers+nOuter)) - 10) / normalizef;
                        constraint.gF{1} = @(x) ...
                            obj.GradDynamicProblem(x, nLayers, nOuter, eye(nLayers), S) / normalizef;

                        % constraint.cF{2} = @(x) -constraint.cF{1}(x) - 10/normalizef;
                        % constraint.gF{2} = @(x) -constraint.gF{1}(x);
                    end
            end

            s.type    = "fmincon";
            s.maxIter = 60;
            s.tolerance = 1e-5;

            switch optimCase
                case 'STATIC',  s.constraintCase = {'INEQUALITY'};
                % case 'DYNAMIC', s.constraintCase = {'INEQUALITY','INEQUALITY'};
                case 'DYNAMIC', s.constraintCase = {'INEQUALITY'};
            end

            cParams.cost         = cost;
            cParams.constraint   = constraint;
            cParams.initialGuess = x0_opt;
            cParams.settings     = s;
            cParams.printingPath = true;
            problem = AcademicProblem(cParams);
            fprintf('Starting optimization (%s)...\n', optimCase);
            tOpt = tic;
            problem.compute();
            tOptElapsed = toc(tOpt);
            fprintf('  Optimization finished in %.2f s\n\n', tOptElapsed);

            xStar = problem.result.fun.fValues;

            if isSymmetric && forceSymmetry
                h_vals           = T * xStar(1:nHalf);
                dTheta_outer_opt = xStar(nHalf+1:nHalf+nOuter);
                deltaTheta_opt   = T * S * dTheta_outer_opt;
            else
                h_vals           = xStar(1:nLayers);
                dTheta_outer_opt = xStar(nLayers+1:nLayers+nOuter);
                deltaTheta_opt   = S * dTheta_outer_opt;
            end

            fprintf('\n==================================================\n');
            fprintf('               OPTIMIZATION RESULTS\n');
            fprintf('==================================================\n');
            for i = 1:nLayers
                fprintf('Layer %2d (%s): h=%8.6f m  |  dTheta=%8.4f deg  |  FinalTheta=%8.4f deg\n', ...
                    i, obj.materialLayers{i}, h_vals(i), deltaTheta_opt(i), ...
                    obj.rotation(i) + deltaTheta_opt(i));
            end
            fprintf('--------------------------------------------------\n');
            fprintf('Total Thickness    : %10.6f m  (Limit: %.4f m)\n', sum(h_vals), h_max);
            fprintf('Initial deltaTheta : %10.6f deg\n', deltaTheta0);
            fprintf('Outer deltaThethas : ');
            fprintf('%8.4f deg  ', dTheta_outer_opt);
            fprintf('\n');

            switch optimCase
                case 'STATIC'
                    finalCost = obj.staticProblem(h_vals, deltaTheta_opt);
                    fprintf('Final Displacement : %10.6e\n', finalCost);
                case 'DYNAMIC'
                    finalCost = obj.dynamicProblem(h_vals, deltaTheta_opt);
                    fprintf('Final frequency    : %10.6e\n', finalCost);
            end

            finalMass = obj.computeMass(h_vals);
            fprintf('Final Mass (kg)    : %10.6e  (Normalized: %10.6e)\n', finalMass, finalMass/m0);
            fprintf('==================================================\n\n');
            % =========================================================================
            % TOTAL ELAPSED TIME
            % =========================================================================
            tTotalElapsed = toc(tTotal);
            [hh, mm, ss]  = obj.formatTime(tTotalElapsed);
            fprintf('--------------------------------------------------\n');
            fprintf('  Optim. time       : %10.2f s\n', tOptElapsed);
            fprintf('  TOTAL run time    : %02d h %02d min %05.2f s\n', hh, mm, ss);
            fprintf('==================================================\n\n');
        end
    end

    methods (Access = private)

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
                    uMesh = UnfittedMesh(sUm);
                    uMesh.compute(ls);
                    obj.mesh = uMesh.createInnerMesh();
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
            g      = GeometricalFunction(gPar);
            phiFun = g.computeLevelSetFunction(mesh);
            ls     = phiFun.fValues;
        end

        function createSolutionField(obj)
            obj.uFun     = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.thetaFun = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.wFun     = LagrangianFunction.create(obj.mesh,1,'P1');
        end

        function wMax = staticProblem(obj, h, deltaTheta)
            obj.createMaterialProperties(h, deltaTheta);
            obj.createBoundaryConditions()
            LHS = obj.createLHS();
            RHS = obj.createRHS();
            x   = LHS\RHS;

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            wF    = x((nU+nTheta+1):(nU+nTheta+nW), 1);
            dofFW = obj.computeFreeDofs(obj.bcW);

            wT = zeros(obj.wFun.nDofs, 1);
            wT(dofFW,1) = wF;
            wT = reshape(wT, [], obj.wFun.ndimf);
            obj.wFun.setFValues(wT);
            wMax = max(wT);
        end

        function freq1 = dynamicProblem(obj, h, deltaTheta)
            obj.createMaterialProperties(h, deltaTheta);
            obj.createBoundaryConditions()
            LHS  = obj.createLHS();
            MLHS = obj.createMassLHS();

            nModes = 3;
            [~, lambda]   = eigs(LHS, MLHS, nModes, 'smallestabs');
            omega_squared = real(diag(lambda));
            omega_squared(omega_squared < 0) = 0;
            freq1 = sqrt(omega_squared(1)) / (2*pi);

            obj.last_h_eval  = h;
            obj.last_f1_eval = freq1;
        end

        function mass = computeMass(obj, h)
            [~,~,~,rho,~] = obj.getMaterialProperties(obj.materialLayers);
            mass = obj.area * sum(rho .* h);
        end

        function g = GradComputeMass(obj)
            [~,~,~,rho,~] = obj.getMaterialProperties(obj.materialLayers);
            g = obj.area .* rho;
        end

        function createMaterialProperties(obj, h, deltaTheta)

            materialName = obj.materialLayers;

            % deltaTheta es un vector nLayers×1 (capas interiores = 0)
            Rotation = obj.rotation + deltaTheta(:);

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
                z0   = z_int(k);   z1 = z_int(k+1);
                C_ps = obj.material.planeStressReduction(C_k);
                Qmatrix(:,:,k) = obj.material.tensorToVoigt(C_ps);

                A_tens = A_tens + C_ps * (z1 - z0);
                B_tens = B_tens + 0.5  * C_ps * (z1^2 - z0^2);
                D_tens = D_tens + 1/3  * C_ps * (z1^3 - z0^3);
                H_tens = H_tens + C_k  * (z1 - z0);
            end

            obj.A_tensor = ConstantFunction.create(A_tens(1:2,1:2,1:2,1:2), obj.mesh);
            obj.B_tensor = ConstantFunction.create(B_tens(1:2,1:2,1:2,1:2), obj.mesh);
            obj.D_tensor = ConstantFunction.create(D_tens(1:2,1:2,1:2,1:2), obj.mesh);

            H_2x2 = [H_tens(1,3,1,3), H_tens(1,3,2,3);
                      H_tens(2,3,1,3), H_tens(2,3,2,3)];
            obj.H_tensor = ConstantFunction.create(H_2x2, obj.mesh);

            index    = [1,2,6];
            A_matrix = obj.material.tensorToVoigt(A_tens);
            B_matrix = obj.material.tensorToVoigt(B_tens);
            D_matrix = obj.material.tensorToVoigt(D_tens);
            F_bar    = A_matrix(4,4) - A_matrix(4,5)^2 / A_matrix(5,5);

            ABD        = [A_matrix(index,index), B_matrix(index,index);
                          B_matrix(index,index), D_matrix(index,index)];
            Compliance = inv(ABD);
            beta1i     = Compliance(1:3,4:6);   beta1i   = beta1i(1,:);
            delta_1i   = Compliance(4:6,4:6);   delta_1i = delta_1i(1,:);

            H_v = zeros(1,nLayers);
            J_v = zeros(1,nLayers);
            for k = 1:nLayers
                H_v(k) = Qmatrix(1,1,k)*beta1i(1)   + Qmatrix(1,2,k)*beta1i(2)   + Qmatrix(1,6,k)*beta1i(3);
                J_v(k) = Qmatrix(1,1,k)*delta_1i(1) + Qmatrix(1,2,k)*delta_1i(2) + Qmatrix(1,6,k)*delta_1i(3);
            end

            sum_term = 0;
            for k = 1:nLayers
                zk = obj.zLayer{k};  zk1 = obj.zLayer{k+1};
                Tk = 0;  Uk = 0;
                for m = 1:(k-1)
                    dz  = obj.zLayer{m+1} - obj.zLayer{m};
                    dz2 = obj.zLayer{m+1}^2 - obj.zLayer{m}^2;
                    Tk  = Tk + H_v(m) * dz;
                    Uk  = Uk + (J_v(m)/2) * dz2;
                end
                Hk = H_v(k);  Jk = J_v(k);

                Pk = Tk^2 + Hk^2*zk^2 - 2*Tk*Hk*zk + Uk^2 + (Jk^2*zk^4)/4 ...
                   - Uk*Jk*zk^2 + 2*Tk*Uk - Tk*Jk*zk^2 ...
                   - 2*Hk*Uk*zk + Hk*Jk*zk^3;
                Rk = 2*Tk*Hk - 2*Hk^2*zk + 2*Hk*Uk - Hk*Jk*zk^2;
                Vk = Hk^2 - (Jk^2*zk^2)/2 + Uk*Jk + Tk*Jk - zk*Hk*Jk;
                Wk = Hk*Jk;
                Xk = Jk^2/4;

                poly_int = Pk*(zk1-zk)       + (Rk/2)*(zk1^2-zk^2) + ...
                           (Vk/3)*(zk1^3-zk^3) + (Wk/4)*(zk1^4-zk^4) + ...
                           (Xk/5)*(zk1^5-zk^5);

                kss      = Qmatrix(4,4,k) - Qmatrix(4,5,k)^2/Qmatrix(5,5,k);
                sum_term = sum_term + poly_int/kss;
            end
            obj.shearCorrectionFactor = (1/F_bar) / sum_term;
        end

        function LHS = createLHS(obj)
            A = obj.A_tensor;
            f = @(u,v) DDP(SymGrad(v), DDP(A, SymGrad(u)));
            Ku = IntegrateLHS(f, obj.uFun, obj.uFun, obj.mesh, 'Domain', 2);
            Ku = reduceMatrix(obj, Ku, obj.bcU, obj.bcU);

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

            nU   = length(obj.computeFreeDofs(obj.bcU));
            nW   = length(obj.computeFreeDofs(obj.bcW));
            beta = obj.shearCorrectionFactor;
            Ztu  = Zut';
            Zuw  = zeros(nU, nW);
            LHS  = [Ku,   Zut,                    Zuw;
                    Ztu,  Ktheta+beta*Mtheta,      beta*Nthetaw;
                    Zuw', beta*Nthetaw',            beta*Kw];
        end

        function RHS = createRHS(obj)
            p = ConstantFunction.create([0 0], obj.mesh);
            m = ConstantFunction.create([0 0], obj.mesh);
            q = ConstantFunction.create(10240,   obj.mesh);

            RHSu     = IntegrateRHS(@(v) DP(p,v), obj.uFun,     obj.mesh, 'Domain', 2);
            RHSu     = obj.reduceVector(RHSu,     obj.bcU);
            RHStheta = IntegrateRHS(@(v) DP(m,v), obj.thetaFun, obj.mesh, 'Domain', 2);
            RHStheta = obj.reduceVector(RHStheta, obj.bcT);
            RHSw     = IntegrateRHS(@(v) q.*v,    obj.wFun,     obj.mesh, 'Domain', 2);
            RHSw     = obj.reduceVector(RHSw,     obj.bcW);

            RHS = [RHSu; RHStheta; RHSw];
        end

        function MLHS = createMassLHS(obj)
            zInterface = cell2mat(obj.zLayer);
            nL   = length(zInterface) - 1;
            Arho = 0;  Brho = 0;  Drho = 0;
            for k = 1:nL
                z0  = zInterface(k);   z1  = zInterface(k+1);
                rho = obj.density.constant(k);
                Arho = Arho + rho*(z1-z0);
                Brho = Brho + 0.5*rho*(z1^2-z0^2);
                Drho = Drho + rho/3*(z1^3-z0^3);
            end

            Mu  = reduceMatrix(obj, IntegrateLHS(@(u,v) Arho*DP(v,u), obj.uFun,     obj.uFun,     obj.mesh,'Domain',2), obj.bcU, obj.bcU);
            Mut = reduceMatrix(obj, IntegrateLHS(@(u,v) Brho*DP(v,u), obj.uFun,     obj.thetaFun, obj.mesh,'Domain',2), obj.bcU, obj.bcT);
            Mt  = reduceMatrix(obj, IntegrateLHS(@(u,v) Drho*DP(v,u), obj.thetaFun, obj.thetaFun, obj.mesh,'Domain',2), obj.bcT, obj.bcT);
            Mw  = reduceMatrix(obj, IntegrateLHS(@(u,v) Arho*v.*u,    obj.wFun,     obj.wFun,     obj.mesh,'Domain',2), obj.bcW, obj.bcW);

            nU = length(obj.computeFreeDofs(obj.bcU));
            nT = length(obj.computeFreeDofs(obj.bcT));
            nW = length(obj.computeFreeDofs(obj.bcW));

            MLHS = [Mu,          Mut,         zeros(nU,nW);
                    Mut.',       Mt,           zeros(nT,nW);
                    zeros(nW,nU), zeros(nW,nT), Mw];
        end

        function createBoundaryConditions(obj)
            obj.bcU = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcT = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcW = obj.createGeneralBoundaryConditions([1]);
        end

        function bc = createGeneralBoundaryConditions(obj, direct)
            TOL  = 1e-12;
            xMin = min(obj.mesh.coord(:,1));
            sDir{1}.domain    = @(coor) abs(coor(:,1)-xMin) < TOL;
            sDir{1}.direction = direct;
            sDir{1}.value     = 0;
            sDir{1}.ndim      = length(direct);
            dirichletFun = [];
            for i = 1:numel(sDir)
                dirichletFun = [dirichletFun, DirichletCondition(obj.mesh, sDir{i})];
            end
            s.dirichletFun = dirichletFun;
            s.periodicFun  = [];
            s.pointloadFun = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

        function fD = computeFreeDofs(obj, bC)
            fD = setdiff(1:bC.dirichletFun.nDofs, bC.dirichlet_dofs);
        end

        function R = reduceMatrix(obj, M, bcV, bcU)
            R = M(obj.computeFreeDofs(bcV), obj.computeFreeDofs(bcU));
        end

        function R = reduceVector(obj, V, bc)
            R = V(obj.computeFreeDofs(bc), 1);
        end

        function grad = GradDynamicProblem(obj, x, nHalf, nOuter, T, S)
            dh     = 1e-6;
            dtheta = 1.0;

            h_half       = x(1:nHalf);
            dTheta_outer = x(nHalf+1:nHalf+nOuter);

            grad_h = zeros(nHalf,1);
            for i = 1:nHalf
                hf = h_half; hb = h_half;
                hf(i) = hf(i) + dh;  hb(i) = hb(i) - dh;
                grad_h(i) = (obj.dynamicProblem(T*hf, T*S*dTheta_outer) - ...
                             obj.dynamicProblem(T*hb, T*S*dTheta_outer)) / (2*dh);
            end

            grad_theta = zeros(nOuter,1);
            for i = 1:nOuter
                tf = dTheta_outer; tb = dTheta_outer;
                tf(i) = tf(i) + dtheta;  tb(i) = tb(i) - dtheta;
                grad_theta(i) = (obj.dynamicProblem(T*h_half, T*S*tf) - ...
                                 obj.dynamicProblem(T*h_half, T*S*tb)) / (2*dtheta);
            end

            grad = [grad_h; grad_theta];
        end

        function grad = GradStaticProblem(obj, x, nHalf, nOuter, T, S)
            dh     = 1e-6;
            dtheta = 1.0;

            h_half       = x(1:nHalf);
            dTheta_outer = x(nHalf+1:nHalf+nOuter);

            grad_h = zeros(nHalf,1);
            for i = 1:nHalf
                hf = h_half; hb = h_half;
                hf(i) = hf(i) + dh;  hb(i) = hb(i) - dh;
                grad_h(i) = (obj.staticProblem(T*hf, T*S*dTheta_outer) - ...
                             obj.staticProblem(T*hb, T*S*dTheta_outer)) / (2*dh);
            end

            grad_theta = zeros(nOuter,1);
            for i = 1:nOuter
                tf = dTheta_outer; tb = dTheta_outer;
                tf(i) = tf(i) + dtheta;  tb(i) = tb(i) - dtheta;
                grad_theta(i) = (obj.staticProblem(T*h_half, T*S*tf) - ...
                                 obj.staticProblem(T*h_half, T*S*tb)) / (2*dtheta);
            end

            grad = [grad_h; grad_theta];
        end

        function [E, nu, G, density, type] = getMaterialProperties(obj, materialName)
            msi = 6.89476e9;

            db.Aluminum = struct('type','ISOTROPIC', 'E',10.6*msi, 'nu',0.33, 'G',3.38*msi,  'density',2700);
            db.Copper   = struct('type','ISOTROPIC', 'E',18.0*msi, 'nu',0.33, 'G',6.39*msi,  'density',8960);
            db.Steel    = struct('type','ISOTROPIC', 'E',30.0*msi, 'nu',0.29, 'G',11.24*msi, 'density',7850);
            db.AS   = struct('type','ORTHOTROPIC','E',[20.0,1.3,1.3]*msi, 'nu',[0.30,0.30,0.49],'G',[1.03,1.03,0.90]*msi,'density',1600);
            db.EpT  = struct('type','ORTHOTROPIC','E',[19.0,1.5,1.5]*msi, 'nu',[0.22,0.22,0.49],'G',[1.00,0.90,0.90]*msi,'density',1600);
            db.T300_914_C = struct('type','ORTHOTROPIC','E',[138.0,11.0,11.0]*1e9, 'nu',[0.28,0.28,0.40],'G',[5.5,5.5,3.928]*1e9,'density',1580);
            db.Ep1  = struct('type','ORTHOTROPIC','E',[7.8,2.6,2.6]*msi,  'nu',[0.25,0.25,0.34],'G',[1.30,1.30,0.50]*msi,'density',1900);
            db.Ep2  = struct('type','ORTHOTROPIC','E',[5.6,1.2,1.3]*msi,  'nu',[0.26,0.26,0.34],'G',[0.60,0.60,0.50]*msi,'density',2000);
            db.BrEp = struct('type','ORTHOTROPIC','E',[30.0,3.0,3.0]*msi, 'nu',[0.30,0.25,0.25],'G',[1.00,1.00,0.60]*msi,'density',2000);

            if ischar(materialName), materialName = {materialName}; end
            nMat = numel(materialName);

            fieldNames = cellfun(@(s) strrep(strrep(s,'-','_'),' ','_'), ...
                materialName, 'UniformOutput', false);

            for i = 1:nMat
                if ~isfield(db, fieldNames{i})
                    error('Material "%s" not found.', materialName{i});
                end
            end

            hasOrtho = any(cellfun(@(f) strcmp(db.(f).type,'ORTHOTROPIC'), fieldNames));
            if hasOrtho, type = 'ORTHOTROPIC'; else, type = 'ISOTROPIC'; end

            if strcmp(type,'ISOTROPIC')
                E = zeros(nMat,1); nu = zeros(nMat,1); G = zeros(nMat,1); density = zeros(nMat,1);
                for i = 1:nMat
                    E(i) = db.(fieldNames{i}).E;  nu(i) = db.(fieldNames{i}).nu;
                    G(i) = db.(fieldNames{i}).G;  density(i) = db.(fieldNames{i}).density;
                end
            else
                E = zeros(nMat,3); nu = zeros(nMat,3); G = zeros(nMat,3); density = zeros(nMat,1);
                for i = 1:nMat
                    m = db.(fieldNames{i});
                    if strcmp(m.type,'ISOTROPIC')
                        E(i,:) = [m.E,m.E,m.E]; nu(i,:) = [m.nu,m.nu,m.nu]; G(i,:) = [m.G,m.G,m.G];
                    else
                        E(i,:) = m.E; nu(i,:) = m.nu; G(i,:) = m.G;
                    end
                    density(i) = m.density;
                end
            end

            if nMat == 1 && size(E,2) > 1
                E = reshape(E,1,[]); nu = reshape(nu,1,[]); 
                G = reshape(G,1,[]); density = reshape(density,1,[]);
            end
        end

        function isSymmetric = checkLaminateSymmetry(obj, Rotation)
            nL = length(obj.materialLayers);
            isSymmetric = true;
            for i = 1:floor(nL/2)
                opp = nL - i + 1;
                if ~strcmp(obj.materialLayers{i}, obj.materialLayers{opp})
                    isSymmetric = false;
                end
                if abs(Rotation(i) - Rotation(opp)) > 1e-6
                    isSymmetric = false;
                end
            end
        end

        % =========================================================================
        % TIME FORMATTER
        % =========================================================================
        function [hh, mm, ss] = formatTime(~, t)
            hh = floor(t / 3600);
            mm = floor(mod(t, 3600) / 60);
            ss = mod(t, 60);
        end

    end
end