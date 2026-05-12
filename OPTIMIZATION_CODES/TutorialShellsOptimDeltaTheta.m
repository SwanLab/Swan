classdef TutorialShellsOptimDeltaTheta < handle

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
        function obj = TutorialShellsOptimDeltaTheta()

            clc; close all;

            tic;
            optimCase = 'DYNAMIC';
            % STATIC
            % DYNAMIC

            % Laminate Inputs 
            obj.materialLayers = {'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT';
                'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'};

            Rotation = [0; 0; 45; -45; 90; 90; -45; 45; 0; 0];  % degrees
            Rotation = 25*ones(size(obj.materialLayers)) + Rotation;
            obj.rotation = Rotation;

            % Symmetry inputs 
            forceSymmetry = true; 

            isSymmetric = obj.checkLaminateSymmetry(obj.rotation);

            % Initialization 
            h_max = 0.5; 
            nLayers = length(obj.materialLayers);
            h0 = h_max * ones(nLayers,1)/nLayers;    
            deltaTheta0 = 15; 

            obj.createMesh('wingShape')  % 'triangleMesh' // 'wingShape'
            obj.createSolutionField()

            % Boundary and normalization values 
            minDeltaTheta = -89.9;
            maxDeltaTheta = 89.9; 
            normalizef = 20; 
            maxval = obj.staticProblem(h0, 0);  
            m0     = obj.computeMass(h0);
            switch optimCase
                case 'STATIC'
                    if isSymmetric && forceSymmetry
                        nHalf = ceil(nLayers/2);
                        x0_opt = [h0(1:nHalf); deltaTheta0];

                        I = eye(nHalf);
                        if mod(nLayers,2) == 0
                            T = [I; flipud(I)];            % Even
                        else
                            T = [I; flipud(I(1:end-1,:))]; % Odd
                        end

                        % Reduced bounds (h + deltaTheta)
                        s.ub = [h_max * ones(nHalf,1);  maxDeltaTheta];
                        s.lb = [0.01  * ones(nHalf,1);   minDeltaTheta];

                        cost.cF = @(x) obj.computeMass(T * x(1:nHalf)) / m0;
                        cost.gF = @(x) [T' * obj.GradComputeMass() / m0; 0];

                        constraint.cF{1} = @(x)  (obj.staticProblem(T * x(1:nHalf), x(nHalf+1)) / maxval - 0.75);
                        constraint.gF{1} = @(x) obj.GradStaticProblem(x, nHalf, T) / maxval;

                    else
                        x0_opt = [h0; deltaTheta0];

                        s.ub = [h_max * ones(nLayers,1);  maxDeltaTheta];
                        s.lb = [0.01  * ones(nLayers,1);   minDeltaTheta];

                        cost.cF = @(x) obj.computeMass(x(1:nLayers)) / m0;
                        cost.gF = @(x) [obj.GradComputeMass() / m0; 0];

                        constraint.cF{1} = @(x)  (obj.staticProblem(x(1:nLayers), x(nLayers+1)) / maxval - 0.75);
                        constraint.gF{1} = @(x) obj.GradStaticProblem(x, nLayers, eye(nLayers)) / maxval;
                    end

                case 'DYNAMIC'
                    if isSymmetric && forceSymmetry
                        nHalf = ceil(nLayers/2);
                        x0_opt = [h0(1:nHalf); deltaTheta0];

                        I = eye(nHalf);
                        if mod(nLayers,2) == 0
                            T = [I; flipud(I)];            % Even
                        else
                            T = [I; flipud(I(1:end-1,:))]; % Odd
                        end

                        % Reduced bounds (h + deltaTheta)
                        s.ub = [h_max * ones(nHalf,1);  maxDeltaTheta];
                        s.lb = [0.01  * ones(nHalf,1);   minDeltaTheta];

                        cost.cF = @(x) obj.computeMass(T * x(1:nHalf)) / m0;
                        cost.gF = @(x) [T' * obj.GradComputeMass() / m0; 0];

                        constraint.cF{1} = @(x)  (obj.dynamicProblem(T * x(1:nHalf), x(nHalf+1)) - 20) / normalizef;
                        constraint.gF{1} = @(x) obj.GradDynamicProblem(x, nHalf, T) / normalizef;

                        constraint.cF{2} = @(x)  -constraint.cF{1}(x) - 10 / normalizef;
                        constraint.gF{2} = @(x)  -constraint.gF{1}(x);

                    else
                        x0_opt = [h0; deltaTheta0];

                        s.ub = [h_max * ones(nLayers,1);  maxDeltaTheta];
                        s.lb = [0.01  * ones(nLayers,1);   minDeltaTheta];

                        cost.cF = @(x) obj.computeMass(x(1:nLayers)) / m0;
                        cost.gF = @(x) [obj.GradComputeMass() / m0; 0]; % masa no depende de deltaTheta

                        constraint.cF{1} = @(x)  (obj.dynamicProblem(x(1:nLayers), x(nLayers+1)) - 20);
                        constraint.gF{1} = @(x) obj.GradDynamicProblem(x, nLayers, eye(nLayers));

                        constraint.cF{2} = @(x)  -constraint.cF{1}(x) - 10 / normalizef;
                        constraint.gF{2} = @(x)  -constraint.gF{1}(x);
                    end
            end

            s.type           = "fmincon";
            s.maxIter        = 17;

            switch optimCase
                case 'STATIC'
                    s.constraintCase = {'INEQUALITY'};
                case 'DYNAMIC'
                    s.constraintCase = {'INEQUALITY','INEQUALITY'};
            end
            cParams.cost         = cost;
            cParams.constraint   = constraint;
            cParams.initialGuess = x0_opt;
            cParams.settings     = s;
            cParams.printingPath = true;
            problem              = AcademicProblem(cParams);
            problem.compute();
            computingTime = toc;

            xStar      = problem.result.fun.fValues;
            
            if isSymmetric && forceSymmetry
                h_vals     = T * xStar(1:nHalf);
                deltaTheta_opt = xStar(nHalf + 1);
            else
                h_vals = xStar(1:nLayers);
                deltaTheta_opt = xStar(end);
            end
            
            fprintf('\n==================================================\n');
            fprintf('               OPTIMIZATION RESULTS\n');
            fprintf('==================================================\n');
            
            % Print the thickness of each layer
            for i = 1:length(h_vals)
                fprintf('Layer %2d (%s): %10.6f m\n', i, obj.materialLayers{i}, h_vals(i));
            end
            
            fprintf('--------------------------------------------------\n');
            fprintf('Total Thickness : %10.6f m  (Limit: %.4f m)\n', sum(h_vals), h_max);
            fprintf('Optimal deltaTheta       : %10.6f deg\n', deltaTheta_opt);
            
            % Calculate and show the final cost/deflection
            switch optimCase
                case 'STATIC'
                    finalCost = obj.staticProblem(h_vals, deltaTheta_opt);
                    fprintf('Final Displacement        : %10.6e\n', finalCost);
                case 'DYNAMIC'
                    finalCost = obj.dynamicProblem(h_vals, deltaTheta_opt);
                    fprintf('Final frequency           : %10.6e\n', finalCost);
            end
            % Print final mass
            finalMass = obj.computeMass(h_vals);
            fprintf('Final Mass (kg)           : %10.6e  (Normalized: %10.6e)\n', finalMass, finalMass / m0);
            fprintf('==================================================\n\n');
            fprintf('Total computing time (min)  : %10.6f\n', computingTime/60);


            
        end
    end

    methods (Access = private)

        function createMesh(obj,meshtype)

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


        function createSolutionField(obj)
           obj.uFun     = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.thetaFun = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.wFun     = LagrangianFunction.create(obj.mesh,1,'P1');
        end
        
        %% staticProblem
        function wMax = staticProblem(obj, h, deltaTheta)

            obj.createMaterialProperties(h,deltaTheta);
            
            obj.createBoundaryConditions()
            LHS = obj.createLHS();
            RHS = obj.createRHS();

            x = LHS\RHS;

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            wF = x((nU+nTheta+1):(nU+nTheta+nW),1);
            dofFW = obj.computeFreeDofs(obj.bcW);

            wT = zeros(obj.wFun.nDofs,1);
            wT(dofFW,1) = wF; 
            wT = reshape(wT,[], obj.wFun.ndimf);
            obj.wFun.setFValues(wT);

            wMax = max(wT);

            % Minimize strain energy 
            % wMax = RHS.'*x;

        end
        
        %% dynamicProblem
        function freq1 = dynamicProblem(obj, h, deltaTheta)

            obj.createMaterialProperties(h,deltaTheta);
            
            obj.createBoundaryConditions()
            LHS = obj.createLHS();
            MLHS = obj.createMassLHS();

            nModes = 3;

            [modes, lambda] = eigs(LHS, MLHS, nModes, 'smallestabs');
            omega_squared = diag(lambda);
            omega_squared = real(omega_squared);
            omega_squared(omega_squared < 0) = 0;  
            omega_n = sqrt(omega_squared);  % rad/s
            
            omega1 = omega_n(1);   
            freq1 = omega1 / 2 / pi;

            obj.last_h_eval = h;
            obj.last_f1_eval = freq1;
        end

        function mass = computeMass(obj,h)
            materialName = obj.materialLayers;
            [~,~,~, rho, ~] = obj.getMaterialProperties(materialName);
            mass = obj.area * sum(rho .* h);
        end

        function gradienteCostFunc = GradComputeMass(obj)
            materialName = obj.materialLayers;
            [~,~,~, rho, ~] = obj.getMaterialProperties(materialName);
            Area = obj.area; 

            gradienteCostFunc = Area .* rho;
        end

        %% createMaterialProperties
        function createMaterialProperties(obj, h, deltaTheta)
            % E = 30;
            % nu = 0.3;
            % Density = 0.1;
            % G = E/ 2*(1+nu);
            % obj.young   = ConstantFunction.create(E, obj.mesh);
            % obj.poisson = ConstantFunction.create(nu, obj.mesh);
            % obj.density = ConstantFunction.create(Density, obj.mesh);
            % obj.area    = ConstantFunction.create(h,obj.mesh);
            % obj.shear   = ConstantFunction.create(G,obj.mesh);
            % obj.inertia = ConstantFunction.create(h^3*1/12,obj.mesh);

            
            materialName = obj.materialLayers;

            Rotation = obj.rotation;
            Tred = [1; 1; 1; 1; 1]; 
            T = [Tred; flipud(Tred)];
            Rotation = T*deltaTheta + Rotation;

            % Get material properties from database
            [E, nu, G, rho, type] = obj.getMaterialProperties(materialName);
            obj.materialProperties = type;
         
            % =========================================================================
            % STORE PROPERTIES IN OBJECT
            % =========================================================================

            obj.young      = ConstantFunction.create(E, obj.mesh);
            obj.poisson    = ConstantFunction.create(nu, obj.mesh);
            obj.shear      = ConstantFunction.create(G, obj.mesh);
            obj.thickness  = ConstantFunction.create(h, obj.mesh);
            obj.inertia    = ConstantFunction.create(1, obj.mesh);  
            obj.density    = ConstantFunction.create(rho,obj.mesh);

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

        function RHS = createRHS(obj)
            p = ConstantFunction.create([0 0],obj.mesh);
            m = ConstantFunction.create([0 0],obj.mesh);
            q = ConstantFunction.create(100,obj.mesh);

            fu = @(v) DP(p,v);
            RHSu = IntegrateRHS(fu,obj.uFun,obj.mesh,'Domain',2);
            RHSu = obj.reduceVector(RHSu,obj.bcU);

            ftheta   = @(v) DP(m,v);
            RHStheta = IntegrateRHS(ftheta,obj.thetaFun,obj.mesh,'Domain',2);
            RHStheta = obj.reduceVector(RHStheta,obj.bcT);

            fw = @(v) q.*v;
            RHSw = IntegrateRHS(fw,obj.wFun,obj.mesh,'Domain',2);
            RHSw = obj.reduceVector(RHSw,obj.bcW);
            

            RHS = [RHSu;RHStheta;RHSw];
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


        function createBoundaryConditions(obj)
            obj.bcU = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcT = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcW = obj.createGeneralBoundaryConditions([1]);
        end

        function bc = createGeneralBoundaryConditions(obj,direct)
            TOL = 1e-12;
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            isLeft  = @(coor)  abs(coor(:,1)-xMin)< TOL;
            %isRight = @(coor)  abs(coor(:,1)-xMax)< TOL;
            isTop   = @(coor)  abs(coor(:,2)-yMax)< TOL;
            isBotom = @(coor)  abs(coor(:,2)-yMin)< TOL;

            sDir{1}.domain    = @(coor) isLeft(coor);
            %sDir{1}.domain    = @(coor) isLeft(coor) | isRight(coor) | isTop(coor) | isBotom(coor);
            
            sDir{1}.direction = direct;
            sDir{1}.value     = 0;                    
            sDir{1}.ndim = length(direct);


            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            s.periodicFun  = [];
            s.pointloadFun    = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);            
        end

        function fD = computeFreeDofs(obj,bC)
            dofs = 1:bC.dirichletFun.nDofs;
            fD   = setdiff(dofs,bC.dirichlet_dofs); 
        end

        function LHSred = reduceMatrix(obj,LHS,bcV,bcU)
            fdofV = obj.computeFreeDofs(bcV);
            fdofU = obj.computeFreeDofs(bcU);
            LHSred = LHS(fdofV,fdofU);
        end

        function RHSred = reduceVector(obj,RHS,bc)
            fdofV = obj.computeFreeDofs(bc);
            RHSred = RHS(fdofV,1);
        end

        %% Gradients
        function grad = GradDynamicProblem(obj, x, nHalf, T)
            dh     = 1e-6;
            dtheta = 1;

            h_half     = x(1:nHalf);
            deltaTheta = x(nHalf + 1);

            grad_h = zeros(nHalf, 1);

            for i = 1:nHalf
                h_f = h_half; h_b = h_half;
                h_f(i) = h_half(i) + dh;
                h_b(i) = h_half(i) - dh;
                f2 = obj.dynamicProblem(T * h_f, deltaTheta);  % expande con T aquí
                f1 = obj.dynamicProblem(T * h_b, deltaTheta);
                grad_h(i) = (f2 - f1) / (2*dh);
            end

            f2 = obj.dynamicProblem(T * h_half, deltaTheta + dtheta);
            f1 = obj.dynamicProblem(T * h_half, deltaTheta - dtheta);

            grad = [grad_h; (f2 - f1) / (2*dtheta)];
        end

        function grad = GradStaticProblem(obj, x, nHalf, T)
            dh     = 1e-6;
            dtheta = 1;

            h_half     = x(1:nHalf);
            deltaTheta = x(nHalf + 1);

            grad_h = zeros(nHalf, 1);

            for i = 1:nHalf
                h_f = h_half; h_b = h_half;
                h_f(i) = h_half(i) + dh;
                h_b(i) = h_half(i) - dh;
                f2 = obj.staticProblem(T * h_f, deltaTheta);
                f1 = obj.staticProblem(T * h_b, deltaTheta);
                grad_h(i) = (f2 - f1) / (2*dh);
            end

            f2 = obj.staticProblem(T * h_half, deltaTheta + dtheta);
            f1 = obj.staticProblem(T * h_half, deltaTheta - dtheta);

            grad = [grad_h; (f2 - f1) / (2*dtheta)];
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


        function isSymmetric = checkLaminateSymmetry(obj, Rotation)
            nLayers = length(obj.materialLayers);
            isSymmetric = true;
            tol = 1e-6; 

            for i = 1:floor(nLayers/2)
                opposite_i = nLayers - i + 1;

                % Check material 
                if ~strcmp(obj.materialLayers{i}, obj.materialLayers{opposite_i})
                    isSymmetric = false;
                end

                % Check orientation 
                if abs(Rotation(i) - Rotation(opposite_i)) > tol
                    isSymmetric = false;
                end

                % % Check thickness
                % if abs(h(i) - h(opposite_i)) > tol
                %     isSymmetric = false;
                % end
            end

            % if isSymmetric
            %     fprintf('Symmetric laminate.\n');
            % else
            %     fprintf('Non-symmetric laminate.\n');
            % end
            % fprintf('-----------------------------------------\n');
        end


    end

end