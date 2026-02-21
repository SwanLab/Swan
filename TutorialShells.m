classdef TutorialShells < handle


    properties (Access=private)
    % properties (Access = {?TutorialShells, ?auxiliaryFunctions})  %(Access = private)
        mesh
        young
        area
        shear
        inertia
        poisson
        uFun
        thetaFun
        wFun
        bcU,bcT,bcW
        lhs,RHSq
        solverType
        type, values, fun
        bcCase
        A_coupling, D_bending, B_coupling, H_shear, zLayer, interfaceIndex
        A_Matrix, B_Matrix, D_Matrix, H_Matrix 
        materialProperties, materialLayers
    end

    methods (Access = public)
        %% TutorialShells
        function obj = TutorialShells()
            clc; close all; 
            
            obj.createMesh()
            obj.createMaterialProperties()
            obj.createSolutionField()

            obj.solverType = 'REDUCED';


            obj.createBoundaryConditions();                    % BOUNDARY CONDITIONS 
            LHS = obj.createLHS();
            obj.lhs = LHS;
            RHS = obj.createRHS();          % BC --> Distributed forces 

            x = LHS\RHS;

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));            
            uF = x(1:nU,1);
            tF = x((nU+1):(nU+nTheta),1);
            wF = x((nU+nTheta+1):(nU+nTheta+nW),1);

            dofFT = obj.computeFreeDofs(obj.bcT);
            dofFU = obj.computeFreeDofs(obj.bcU);
            dofFW = obj.computeFreeDofs(obj.bcW);

            uT = zeros(obj.uFun.nDofs,1);
            uT(dofFU,1) = uF; 
            uT = reshape(uT,obj.uFun.ndimf,[])';
            obj.uFun.setFValues(uT);

            wT = zeros(obj.wFun.nDofs,1);
            wT(dofFW,1) = wF; 
            wT = reshape(wT,[], obj.wFun.ndimf);
            obj.wFun.setFValues(wT);
            
            thetaT = zeros(obj.thetaFun.nDofs,1);    
            thetaT(dofFT,1) = tF; 
            thetaT = reshape(thetaT,obj.thetaFun.ndimf,[])';
            obj.thetaFun.setFValues(thetaT);

            h = obj.zLayer;
            
            % Layer medium-plane 
            for i = 1:numel(h)-1
                z{i} = 0.5*(h{i+1}+h{i});
            end

            % Obtain epsilon U, Theta, transverse strains 
            [epsilonU_nodal, epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal] = obj.createEpsilons();

            % Obtain Strain and Stress Functions on a certain height 
            [strainFun, stressFun] = obj.createStrainStressFunctions(h, epsilonU_nodal, epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal);

            % Obtain Internal Forces 
            [Nfun, Mfun, Qfun] = obj.internalForces(epsilonU_nodal,epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal);


            %% Save the solution fValues of uFun, thetaFun and wFun to MAT-files
            ufun_vals = obj.uFun.fValues;        % nodal DOF values for u (matrix size: nNodes x ndimf)
            theta_vals = obj.thetaFun.fValues;   % nodal DOF values for theta (matrix size: nNodes x ndimf)
            w_vals = obj.wFun.fValues;           % nodal DOF values for w (matrix size: nNodes x ndimf)

            % % Use full paths or simple filenames in current folder
            % save('uFun_fValues.mat','ufun_vals');
            % save('thetaFun_fValues.mat','theta_vals');
            % save('wFun_fValues.mat','w_vals');



            % Load saved solution files if they exist and compare to current solutions
            savedDiffs = struct();

            % Compare uFun
            if exist('uFun_fValues.mat','file') == 2
                S = load('uFun_fValues.mat','ufun_vals');
                if isfield(S,'ufun_vals')
                    diffU = obj.uFun.fValues - S.ufun_vals;
                    savedDiffs.u = diffU;

                    norm_diff = norm(diffU(:));
                    norm_u = norm(obj.uFun.fValues(:));
                    rel_error = norm_diff / norm_u;

                    fprintf('=== uFun Comparison ===\n');
                    fprintf('Absolute error (L2 norm): %.2e\n', norm_diff);
                    fprintf('Relative error: %.2e\n', rel_error);

                    if rel_error < 1e-10
                        fprintf('✓ Solutions are identical (numerical noise only)\n');
                    elseif rel_error < 1e-6
                        fprintf('✓ Solutions are equivalent (acceptable numerical differences)\n');
                    elseif rel_error < 1e-3
                        fprintf('⚠ Solutions are similar (minor differences)\n');
                    else
                        fprintf('⚠ Solutions differ significantly\n');
                    end

                    idxU = find(abs(diffU) > 1e-6);
                    if ~isempty(idxU)
                        fprintf('Nonzero differences (>1e-6): %d out of %d\n', numel(idxU), numel(diffU));
                        fprintf('Max absolute difference: %.2e\n', max(abs(diffU(:))));
                    end
                    fprintf('\n');
                end
            else
                fprintf('File uFun_fValues.mat not found.\n');
            end

            % Compare thetaFun
            if exist('thetaFun_fValues.mat','file') == 2
                S = load('thetaFun_fValues.mat','theta_vals');
                if isfield(S,'theta_vals')
                    diffTheta = obj.thetaFun.fValues - S.theta_vals;
                    savedDiffs.theta = diffTheta;

                    norm_diff = norm(diffTheta(:));
                    norm_theta = norm(obj.thetaFun.fValues(:));
                    rel_error = norm_diff / norm_theta;

                    fprintf('=== thetaFun Comparison ===\n');
                    fprintf('Absolute error (L2 norm): %.2e\n', norm_diff);
                    fprintf('Relative error: %.2e\n', rel_error);

                    if rel_error < 1e-10
                        fprintf('✓ Solutions are identical (numerical noise only)\n');
                    elseif rel_error < 1e-6
                        fprintf('✓ Solutions are equivalent (acceptable numerical differences)\n');
                    elseif rel_error < 1e-3
                        fprintf('⚠ Solutions are similar (minor differences)\n');
                    else
                        fprintf('⚠ Solutions differ significantly\n');
                    end

                    idxT = find(abs(diffTheta) > 1e-6);
                    if ~isempty(idxT)
                        fprintf('Nonzero differences (>1e-6): %d out of %d\n', numel(idxT), numel(diffTheta));
                        fprintf('Max absolute difference: %.2e\n', max(abs(diffTheta(:))));
                    end
                    fprintf('\n');
                end
            else
                fprintf('File thetaFun_fValues.mat not found.\n');
            end

            % Compare wFun
            if exist('wFun_fValues.mat','file') == 2
                S = load('wFun_fValues.mat','w_vals');
                if isfield(S,'w_vals')
                    diffW = obj.wFun.fValues - S.w_vals;
                    savedDiffs.w = diffW;

                    norm_diff = norm(diffW(:));
                    norm_w = norm(obj.wFun.fValues(:));
                    rel_error = norm_diff / norm_w;

                    fprintf('=== wFun Comparison ===\n');
                    fprintf('Absolute error (L2 norm): %.2e\n', norm_diff);
                    fprintf('Relative error: %.2e\n', rel_error);

                    if rel_error < 1e-10
                        fprintf('✓ Solutions are identical (numerical noise only)\n');
                    elseif rel_error < 1e-6
                        fprintf('✓ Solutions are equivalent (acceptable numerical differences)\n');
                    elseif rel_error < 1e-3
                        fprintf('⚠ Solutions are similar (minor differences)\n');
                    else
                        fprintf('⚠ Solutions differ significantly\n');
                    end

                    idxW = find(abs(diffW) > 1e-6);
                    if ~isempty(idxW)
                        fprintf('Nonzero differences (>1e-6): %d out of %d\n', numel(idxW), numel(diffW));
                        fprintf('Max absolute difference: %.2e\n', max(abs(diffW(:))));
                    end
                    fprintf('\n');
                end
            else
                fprintf('File wFun_fValues.mat not found.\n');
            end

            
            % Mostrar ordenes de magnitud de las diferencias no nulas guardadas en savedDiffs
            if exist('savedDiffs','var')
                fields = fieldnames(savedDiffs);
                for i=1:numel(fields)
                    name = fields{i};
                    dif = savedDiffs.(name);
                    if isempty(dif)
                        fprintf('No differences recorded for %s.\n', name);
                        continue;
                    end
                    absdif = abs(dif(:));
                    nz = absdif(absdif>0);
                    if isempty(nz)
                        fprintf('All differences for %s are exactly zero.\n', name);
                        continue;
                    end
                    % Compute orders of magnitude (base 10) for nonzero entries
                    orders = floor(log10(nz));
                    unique_orders = sort(unique(orders),'descend');
                    fprintf('Orders of magnitude for nonzero differences in %s:\n', name);
                    for k = 1:numel(unique_orders)
                        o = unique_orders(k);
                        cnt = sum(orders==o);
                        fprintf('  10^{%d}: %d entries (range: %.2e to %.2e)\n', o, cnt, min(nz(orders==o)), max(nz(orders==o)));
                    end
                    % Also print overall stats
                    fprintf('  Total nonzero: %d / %d, min: %.2e, median: %.2e, max: %.2e\n\n', ...
                        numel(nz), numel(absdif), min(nz), median(nz), max(nz));
                end
            else
                fprintf('No saved differences (savedDiffs) available to analyze.\n');
            end

            

            %% PLOT AND PRINT 

            plotMatlab = false;
            printParaview = false;

            % KAPPA ==========================================
            kappa = 1;

            if plotMatlab == true 
                % Displacements
                obj.customPlot(obj.uFun, {'u_{x}', 'u_{y}'});

                % Transverse displacements
                obj.customPlot(obj.wFun, {'w'});

                % Rotations 
                obj.customPlot(obj.thetaFun, {'\theta_{x}', '\theta_{y}'});

                % Strains
                obj.customPlot(strainFun{kappa}, {'\epsilon_{xx}', '\epsilon_{yy}', '\epsilon_{zz}', ...
                    '\epsilon_{xy}', '\epsilon_{xz}', '\epsilon_{yz}'});

                % Stresses
                obj.customPlot(stressFun{1,kappa}, {'\sigma_{xx}', '\sigma_{yy}', '\sigma_{zz}', ...
                    '\sigma_{xy}', '\sigma_{xz}', '\sigma_{yz}'});

                % Resulting forcees 
                obj.customPlot(Nfun, {'N_{xx}', 'N_{yy}', 'N_{xy}'});
                obj.customPlot(Mfun, {'M_{xx}', 'M_{yy}', 'M_{xy}'});
                obj.customPlot(Qfun, {'Q_{xz}', 'Q_{yz}'});
            end

            


            if printParaview == true
                % obj.wFun.print('wfun print','Paraview')
                % obj.uFun.print('ufun print','Paraview')
                % obj.thetaFun.print('thetafun print','Paraview')
                % strainFun{kappa}.print('strain', 'Paraview'); % Kappa defined on plots
                % stressFun{kappa}.print('stress', 'Paraview');
                
                stressFun{1}.print('BOTTOM Layer 1','Paraview')
                stressFun{2}.print('TOP Layer 1','Paraview')
                stressFun{3}.print('BOTTOM Layer 2','Paraview')
                stressFun{4}.print('TOP Layer 2','Paraview')


            end
        end

    end

    % methods (Access = {?TutorialShells, ?auxiliaryFunctions})  % (Access = private)
    methods (Access = private)
        %% createMesh
        function createMesh(obj)
          obj.mesh = UnitTriangleMesh(50,50);
        end

        %% createSolutionField
        function createSolutionField(obj)
           obj.uFun     = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.thetaFun = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.wFun     = LagrangianFunction.create(obj.mesh,1,'P1');
        end

        %% createMaterialProperties
        function createMaterialProperties(obj)
            
            % Material properties [thicknes x direction] direction = 1 for
            % isotropic 
            E = [10; 3; 10];
            nu = [0.3; 0.25; 0.3]; 
            h = [0.25; 0.5; 0.25];
            G = E ./ 2./ (1+nu);
            
            % AS
            % E = [20 1.3 1.3]*1e6*6894.76;
            % nu = [0.3 0.3 0.49];                     % nu_12, nu_13, nu_23
            % h = 1;
            % G = [1.03 1.03 0.9]*1e6*6894.76;         % G_12, G_13, G_23

            obj.materialLayers = 'SINGLE';  % SINGLE / MULTI
            obj.materialProperties = 'ISOTROPIC'; % ISOTROPIC / ANISOTROPIC

            % Set-up 
            obj.young = ConstantFunction.create(E,obj.mesh);
            obj.area = ConstantFunction.create(h,obj.mesh);  % CONSIDERED AS THICKNESS
            obj.shear = ConstantFunction.create(G,obj.mesh);
            obj.inertia = ConstantFunction.create(1,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);

            % Layer properties 
            obj.A_coupling = ConstantFunction.create(sum(E.*h),obj.mesh);
            obj.H_shear = ConstantFunction.create(sum(G.*h),obj.mesh);
            
            D_bend = 0; B_coup = 0;  
            height{1} = -sum(obj.area.constant)/2;
            for i = 2:(length(h)+1)
                height{i} = height{i-1} + h(i-1); 
                D_bend = D_bend + E(i-1)*(height{i}^3-height{i-1}^3)/3;
                B_coup = B_coup + 0.5*E(i-1)*(height{i}^2-height{i-1}^2);
            end
            obj.zLayer = height;
            
            obj.D_bending = ConstantFunction.create(D_bend,obj.mesh);
            obj.B_coupling = ConstantFunction.create(B_coup,obj.mesh); 

            % Matrix layer properties 
            zInterfaces = [obj.zLayer{:}];
            A_matrix = zeros(3);
            B_matrix = zeros(3);
            D_matrix = zeros(3);
            H_matrix = zeros(2);

            % Layer medium-plane (z) as a numeric vector
            z = 0.5 * (zInterfaces(2:end) + zInterfaces(1:end-1));

            for kappa = 1:length(obj.area.constant(:,1))

                ConstitutiveMatrix = obj.createConstitutiveMatrix(z(kappa),zInterfaces,false);
                A_matrix = A_matrix + ConstitutiveMatrix([1,2,4],[1,2,4])*obj.area.constant(kappa);
                B_matrix = B_matrix + 0.5*ConstitutiveMatrix([1,2,4],[1,2,4])*(obj.zLayer{kappa+1}^2-obj.zLayer{kappa}^2);
                D_matrix = D_matrix + 1/3*ConstitutiveMatrix([1,2,4],[1,2,4])*(obj.zLayer{kappa+1}^3-obj.zLayer{kappa}^3);
                H_matrix = H_matrix + ConstitutiveMatrix(5:6,5:6).*obj.area.constant(kappa);

            end

            obj.A_Matrix = A_matrix;
            obj.B_Matrix = B_matrix;
            obj.D_Matrix = D_matrix;
            obj.H_Matrix = H_matrix;        
            
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
            obj.bcCase = 2;     % 1 --> Modified bc // 2 --> Original bc (change q) 
            
            % APPLYED FORCE
            applyedForce = 2;

            % 1 --> Single node
            % 2 --> Right edge

            switch obj.bcCase
                case 1

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

                    

                    switch applyedForce
                        case 1

                            % Apply point load to a single node on the right boundary
                            % Find right boundary nodes
                            rightNodes = find(isRight(obj.mesh.coord));
                            if isempty(rightNodes)
                                error('No nodes found on the right boundary.');
                            end
                            % Choose one node: here pick the middle one (can change as needed)
                            idx = ceil(numel(rightNodes)/2);
                            singleNode = rightNodes(idx);

                            % Define point load structure for that single node
                            sPL{1}.domain    = @(coor) (1:size(coor,1))'==singleNode;
                            sPL{1}.direction = 2;
                            sPL{1}.value     = 1;

                        case 2

                            % Load on Right edge
                            sPL{1}.domain    = @(coor) isRight(coor);
                            sPL{1}.direction = 2;
                            sPL{1}.value     = 1;
                    end

                    pointloadFun = [];
                    for i = 1:numel(sPL)
                        pl = TractionLoad2(obj.mesh, sPL{i}, 'DIRAC');
                        pointloadFun = [pointloadFun, pl];

                    end
                    s.pointloadFun = pointloadFun;


                    s.periodicFun  = [];
                    s.mesh = obj.mesh;
                    bc = BoundaryConditions(s);

                case 2

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
                    s.periodicFun  = [];
                    s.mesh = obj.mesh;
                    bc = BoundaryConditions(s);

            end

        end

        %% createRHS
        function RHS = createRHS(obj)

            p = ConstantFunction.create([0 0],obj.mesh);
            m = ConstantFunction.create([0 0],obj.mesh);
            switch obj.bcCase
                case 1
                    q = ConstantFunction.create(0,obj.mesh);
                case 2
                    q = ConstantFunction.create(1,obj.mesh);
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

            if obj.bcCase == 1
                obj.computeForces();
                RHSq = obj.RHSq;
                RHSq = obj.reduceVector(RHSq,obj.bcW);
                RHSw = RHSw + RHSq;
            end            

            RHS = [RHSu;RHStheta;RHSw];
        end


        %% createLHS
        function LHS = createLHS(obj)
            % A2 = obj.A_coupling;
            % f2 = @(u,v) A2.*DDP(SymGrad(v),SymGrad(u));
            % Ku2 = IntegrateLHS(f2,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            % Ku2 = reduceMatrix(obj,Ku2,obj.bcU,obj.bcU);

            A = obj.A_Matrix; 
            f = @(u,v) TutorialShells.customDDP(A, SymGrad(v), SymGrad(u));            

            Ku = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Ku = reduceMatrix(obj,Ku,obj.bcU,obj.bcU);

            
            % % Subtract Ku2 from Ku and display nonzero entries
            % Kdiff = Ku - Ku2;
            % 
            % % Find nonzero elements (with tolerance)
            % tol = 1e-12;
            % [i, j, s] = find(abs(Kdiff) > tol);
            % 
            % if isempty(s)
            %     disp('No nonzero entries in Ku - Ku2 (within tolerance).');
            % else
            %     fprintf('Found %d nonzero entries in Ku - Ku2\n', length(s));
            %     fprintf('Max absolute difference: %.2e\n', max(abs(s)));
            %     fprintf('First 10 nonzero entries (row, col, value):\n');
            %     for idx = 1:min(10, length(s))
            %         fprintf('%d\t%d\t%.12g\n', i(idx), j(idx), s(idx));
            %     end
            % 
            %     % Calcular error relativo
            %     norm_Kdiff = norm(Kdiff, 'fro');
            %     norm_Ku = norm(Ku, 'fro');
            %     rel_error = norm_Kdiff / norm_Ku;
            % 
            %     fprintf('\nRelative error (Frobenius norm): %.2e\n', rel_error);
            %     if rel_error < 1e-10
            %         fprintf('✓ Matrices are identical\n');
            %     elseif rel_error < 1e-6
            %         fprintf('✓ Matrices are equivalent\n');
            %     else
            %         fprintf('⚠ Matrices differ significantly\n');
            %     end
            % end





            % D = obj.D_bending;
            % f = @(u,v) D.*DDP(SymGrad(v),SymGrad(u));
            D = obj.D_Matrix;
            f = @(u,v) TutorialShells.customDDP(D, SymGrad(v), SymGrad(u));

            Ktheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Ktheta = obj.reduceMatrix(Ktheta,obj.bcT,obj.bcT);

            % H = obj.H_shear; 
            % f = @(u,v) H.*DP(v,u);
            H = obj.H_Matrix;
            f = @(u,v) TutorialShells.customDPvector(H, v, u);

            Mtheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Mtheta = obj.reduceMatrix(Mtheta,obj.bcT,obj.bcT);

            % f = @(u,v) H.*DP(v,Grad(u));
            f = @(u,v) TutorialShells.customDP(H, v, Grad(u));

            Nthetaw = IntegrateLHS(f,obj.thetaFun,obj.wFun,obj.mesh,'Domain',2);            
            Nthetaw = obj.reduceMatrix(Nthetaw,obj.bcT,obj.bcW);

            % f = @(u,v) H.*DP(Grad(v),Grad(u));
            f = @(u,v) TutorialShells.customDPgrad(H, Grad(v), Grad(u));
            
            Kw = IntegrateLHS(f,obj.wFun,obj.wFun,obj.mesh,'Domain',2);  
            Kw = obj.reduceMatrix(Kw,obj.bcW,obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            beta = 1; 

            Zut = zeros(nU,nTheta);
            Zuw = zeros(nU,nW);
            LHS = [Ku Zut Zuw; Zut' (Ktheta+beta*Mtheta) beta*Nthetaw; Zuw' beta*Nthetaw' beta*Kw];
        end

        %% createStrainStressFunctions
        function [strainFun, stressFun] = createStrainStressFunctions(obj, z, epsilonU_nodal, epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal)

            nNodes = size(obj.mesh.coord, 1);
            zInterfaces = [obj.zLayer{:}];

            % Transverse shear strains at nodes
            epsilon_xz_nodal = 0.5 * (dw_dx_nodal + obj.thetaFun.fValues(:,1));
            epsilon_yz_nodal = 0.5 * (dw_dy_nodal + obj.thetaFun.fValues(:,2));

            strainFun = cell(1,numel(z));

            % Check if any z corresponds to an internal interface. 
            thereIsInterface = false;
            interfacePositions = [];
            count = 1; 
            if numel(zInterfaces) > 2
                internalInterfaces = zInterfaces(2:end-1);
                tol = 1e-10;
                for k = 1:numel(z)
                    if any(abs(internalInterfaces - z{k}) <= tol)
                        thereIsInterface = true;
                        interfacePositions(end+1) = k; %#ok<AGROW>
                        fprintf('Found interface at z index: %d\n', k);
                        obj.interfaceIndex(count) = k; 
                        count = count + 1; 
                    end
                end
            end

            if thereIsInterface == false 
                stressTemp = cell(1,numel(z));
                fprintf('Did not found interface\n');
            else
                stressTemp = cell(2,numel(z));
            end

            for i = 1:numel(z)
                strain_nodal = zeros(nNodes, 6);
                strain_nodal(:, 1) = epsilonU_nodal(1, :) + z{i} * epsilonTheta_nodal(1, :);  % exx
                strain_nodal(:, 2) = epsilonU_nodal(2, :) + z{i} * epsilonTheta_nodal(2, :);  % eyy
                strain_nodal(:, 3) = 0;  % ezz = 0
                strain_nodal(:, 4) = epsilonU_nodal(3, :) + z{i} * epsilonTheta_nodal(3, :);  % exy
                strain_nodal(:, 5) = epsilon_xz_nodal;
                strain_nodal(:, 6) = epsilon_yz_nodal;

                % Create LagrangianFunction for strains [exx, eyy, ezz, exy, exz, eyz]
                strainFun{i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
                strainFun{i}.setFValues(strain_nodal);

                % Obtain constitutive matrix for the layer at z{i}
                % If there are multiple layers, check whether current z{i} matches any internal interface (excluding boundaries)
                isInternalInterface = false;
                if numel(zInterfaces) > 2
                    internalInterfaces = zInterfaces(2:end-1);
                    % z{i} may be numeric scalar; use isequal or ismember with tolerance for floating point
                    tol = 1e-10;
                    if any(abs(internalInterfaces - z{i}) <= tol)
                        isInternalInterface = true;
                    end
                end

                C = obj.createConstitutiveMatrix(z{i}, zInterfaces,isInternalInterface);

                if isInternalInterface == false
                    % Obtain stress: sigma = C * epsilon at each node
                    stress_nodal = (C * strain_nodal.').';  % (nNodes x 6)

                    % Create LagrangianFunction for stresses [sxx, syy, szz, sxy, sxz, syz]
                    stressTemp{1,i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
                    stressTemp{1,i}.setFValues(stress_nodal);
                else 
                    for j = 1:2
                        stress_nodal = (C(:,:,j) * strain_nodal.').';  % (nNodes x 6)
                        stressTemp{j,i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
                        stressTemp{j,i}.setFValues(stress_nodal);
                    end
                end
            end


            count = 0;
            fprintf('\n');
            for i = 1:numel(stressTemp(1,:))
                count = count + 1;
                for j = 1:numel(stressTemp(:,1))
                    if ~isempty(stressTemp{j,i})
                        if ~isempty(obj.interfaceIndex) && any(i == obj.interfaceIndex)
                            if j == 1
                                kLayer = find(z{i} >= zInterfaces(1:(end-1)) & z{i} <= zInterfaces(2:end), 1);
                                disp(['Position {' num2str(count) '} of stressFun corresponds to TOP values of layer ' num2str(kLayer)]);
                            else
                                count = count + 1;
                                kLayer = find(z{i} >= zInterfaces(1:(end-1)) & z{i} <= zInterfaces(2:end), 1)+1;
                                disp(['Position {' num2str(count) '} of stressFun corresponds to BOTTOM values of layer ' num2str(kLayer)]);
                            end
                        end
                        stressFun{1,count} = stressTemp{j,i};
                    end
                end
            end
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

        %% createConstitutiveMatrix
        function C = createConstitutiveMatrix(obj,z,z_interfaces,isInterface)
            switch obj.materialProperties
                case 'ISOTROPIC'
                    kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);
                    % when there is layer superposition, takes the value of the bottom layer

                    if isInterface == false
                        nu = obj.poisson.constant(kLayer);
                        lambda = obj.young.constant(kLayer) / (1+nu) / (1-2*nu);
                        % C = lambda* [1-nu, nu, nu, 0, 0, 0;
                        %     nu, 1-nu, nu, 0, 0, 0;
                        %     nu, nu, 1-nu, 0, 0, 0;
                        %     0, 0, 0, 0.5*(1-2*nu), 0, 0;
                        %     0, 0, 0, 0, 0.5*(1-2*nu), 0;
                        %     0, 0, 0, 0, 0, 0.5*(1-2*nu)];
                         C = lambda* [1-nu, nu, nu, 0, 0, 0;
                            nu, 1-nu, nu, 0, 0, 0;
                            nu, nu, 1-nu, 0, 0, 0;
                            0, 0, 0, (1-2*nu), 0, 0;
                            0, 0, 0, 0, (1-2*nu), 0;
                            0, 0, 0, 0, 0, (1-2*nu)];
                    else 
                        for i = 1:2
                            nu = obj.poisson.constant(kLayer);
                            lambda = obj.young.constant(kLayer) / (1+nu) / (1-2*nu);
                            % C(:,:,i) = lambda* [1-nu, nu, nu, 0, 0, 0;
                            %     nu, 1-nu, nu, 0, 0, 0;
                            %     nu, nu, 1-nu, 0, 0, 0;
                            %     0, 0, 0, 0.5*(1-2*nu), 0, 0;
                            %     0, 0, 0, 0, 0.5*(1-2*nu), 0;
                            %     0, 0, 0, 0, 0, 0.5*(1-2*nu)];
                            C(:,:,i) = lambda* [1-nu, nu, nu, 0, 0, 0;
                                nu, 1-nu, nu, 0, 0, 0;
                                nu, nu, 1-nu, 0, 0, 0;
                                0, 0, 0, (1-2*nu), 0, 0;
                                0, 0, 0, 0, (1-2*nu), 0;
                                0, 0, 0, 0, 0, (1-2*nu)];

                            kLayer = kLayer + 1; 
                        end
                    end


                case 'ANISOTROPIC'
                    kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);

                    if isInterface == false
                        % Properties
                        E1 = obj.young.constant(kLayer,1);       E2 = obj.young.constant(kLayer,2);       E3 = obj.young.constant(kLayer,3);
                        nu12 = obj.poisson.constant(kLayer,1);   nu13 = obj.poisson.constant(kLayer,2);   nu23 = obj.poisson.constant(kLayer,3);
                        G12 = obj.shear.constant(kLayer,1);      G13 = obj.shear.constant(kLayer,2);      G23 = obj.shear.constant(kLayer,3);

                        nu21 = nu12*E2/E1;
                        nu32 = nu23*E3/E2;
                        nu31 = nu13*E3/E1;

                        delta = (1-nu12*nu21-nu23*nu32-nu31*nu13-2*nu21*nu32*nu13) / (E1*E2*E3);

                        C11 = (1 - nu23*nu32)/(E2*E3*delta);
                        C12 = (nu21 + nu31*nu23)/(E2*E3*delta);
                        C13 = (nu31 + nu21*nu32)/(E2*E3*delta);
                        C22 = (1 - nu13*nu31)/(E1*E3*delta);
                        C23 = (nu32 + nu12*nu31)/(E1*E3*delta);
                        C33 = (1 - nu12*nu21)/(E1*E2*delta);
                        C44 = G12;
                        C55 = G13;
                        C66 = G23;

                        % C = [ C11 C12 C13 0   0   0;
                        %     C12 C22 C23 0   0   0;
                        %     C13 C23 C33 0   0   0;
                        %     0   0   0   C44 0   0;
                        %     0   0   0   0   C55 0;
                        %     0   0   0   0   0   C66 ];
                        C = [ C11 C12 C13 0   0   0;
                            C12 C22 C23 0   0   0;
                            C13 C23 C33 0   0   0;
                            0   0   0   2*C44 0   0;
                            0   0   0   0   2*C55 0;
                            0   0   0   0   0   2*C66 ];
                    else
                        for i=1:2
                            E1 = obj.young.constant(kLayer,1);       E2 = obj.young.constant(kLayer,2);       E3 = obj.young.constant(kLayer,3);
                            nu12 = obj.poisson.constant(kLayer,1);   nu13 = obj.poisson.constant(kLayer,2);   nu23 = obj.poisson.constant(kLayer,3);
                            G12 = obj.shear.constant(kLayer,1);      G13 = obj.shear.constant(kLayer,2);      G23 = obj.shear.constant(kLayer,3);

                            nu21 = nu12*E2/E1;
                            nu32 = nu23*E3/E2;
                            nu31 = nu13*E3/E1;

                            delta = (1-nu12*nu21-nu23*nu32-nu31*nu13-2*nu21*nu32*nu13) / (E1*E2*E3);

                            C11 = (1 - nu23*nu32)/(E2*E3*delta);
                            C12 = (nu21 + nu31*nu23)/(E2*E3*delta);
                            C13 = (nu31 + nu21*nu32)/(E2*E3*delta);
                            C22 = (1 - nu13*nu31)/(E1*E3*delta);
                            C23 = (nu32 + nu12*nu31)/(E1*E3*delta);
                            C33 = (1 - nu12*nu21)/(E1*E2*delta);
                            C44 = G12;
                            C55 = G13;
                            C66 = G23;

                            C(:,:,i) = [ C11 C12 C13 0   0   0;
                                C12 C22 C23 0   0   0;
                                C13 C23 C33 0   0   0;
                                0   0   0   2*C44 0   0;
                                0   0   0   0   2*C55 0;
                                0   0   0   0   0   2*C66 ];
                            
                            kLayer = kLayer+1;
                        end
                    end
           end
        end

        %% createEpsilons
        function [epsilonU_nodal, epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal] = createEpsilons(obj)
            nNodes = size(obj.mesh.coord, 1);

            % Compute SymGrad of u and theta outside the loop (independent of z)
            strainU     = SymGrad(obj.uFun);
            strainTheta = SymGrad(obj.thetaFun);

            strainUvalues     = strainU.evaluate(obj.mesh);
            strainThetaValues = strainTheta.evaluate(obj.mesh);

            nElem  = size(obj.mesh.connec, 1);

            % Compute gradW at Gauss points and project to nodes
            gradW       = Grad(obj.wFun);
            gradWvalues = gradW.evaluate(obj.mesh);
            dw_dx_gauss = squeeze(gradWvalues(1,:,:));  % (nGauss x nElem)
            dw_dy_gauss = squeeze(gradWvalues(2,:,:));  % (nGauss x nElem)

            % Project gradW to nodes
            dw_dx_nodal = zeros(nNodes, 1);
            dw_dy_nodal = zeros(nNodes, 1);
            node_count  = zeros(nNodes, 1);

            % Project epsilonU to nodes
            epsilonU_nodal     = zeros(nNodes, 3);  % [exx, eyy, exy]
            epsilonTheta_nodal = zeros(nNodes, 3);

            for e = 1:nElem
                elemNodes = obj.mesh.connec(e, :);
                for n = 1:length(elemNodes)
                    nodeIdx = elemNodes(n);
                    dw_dx_nodal(nodeIdx) = dw_dx_nodal(nodeIdx) + mean(dw_dx_gauss(e));
                    dw_dy_nodal(nodeIdx) = dw_dy_nodal(nodeIdx) + mean(dw_dy_gauss(e));

                    epsilonU_nodal(nodeIdx, 1) = epsilonU_nodal(nodeIdx, 1) + mean(squeeze(strainUvalues(1,1,:,e)));
                    epsilonU_nodal(nodeIdx, 2) = epsilonU_nodal(nodeIdx, 2) + mean(squeeze(strainUvalues(2,2,:,e)));
                    epsilonU_nodal(nodeIdx, 3) = epsilonU_nodal(nodeIdx, 3) + mean(squeeze(strainUvalues(1,2,:,e)));

                    epsilonTheta_nodal(nodeIdx, 1) = epsilonTheta_nodal(nodeIdx, 1) + mean(squeeze(strainThetaValues(1,1,:,e)));
                    epsilonTheta_nodal(nodeIdx, 2) = epsilonTheta_nodal(nodeIdx, 2) + mean(squeeze(strainThetaValues(2,2,:,e)));
                    epsilonTheta_nodal(nodeIdx, 3) = epsilonTheta_nodal(nodeIdx, 3) + mean(squeeze(strainThetaValues(1,2,:,e)));

                    node_count(nodeIdx)  = node_count(nodeIdx) + 1;
                end
            end

            dw_dx_nodal = dw_dx_nodal ./ node_count;
            dw_dy_nodal = dw_dy_nodal ./ node_count;

            epsilonU_nodal     = (epsilonU_nodal     ./ node_count).';
            epsilonTheta_nodal = (epsilonTheta_nodal ./ node_count).';
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
        function customPlot(obj, fun, titles)
            % Inputs:
            %   fun    - LagrangianFunction to plot
            %   titles - Cell array of strings with titles for each subplot
            %
            % Example:
            %   customPlot(obj.uFun, {'u_{x}', 'u_{y}'});

            plot(fun);
            ax = findall(gcf, 'type', 'axes');
            ax = flipud(ax);

            for j = 1:length(ax)
                if j <= length(titles)
                    title(ax(j), titles{j});
                end
            end
        end

       
    end

    % En auxiliaryFunctions.m

    methods (Static, Access = public)
        %% customDDP
        function df = customDDP(A, epsV, epsU)

            s.operation = @(xV) TutorialShells.evaluateVoigtProduct(A, epsV, epsU, xV);
            s.mesh = epsV.mesh;
            s.ndimf = 1;
            df = DomainFunction(s);
        end

        function result = evaluateVoigtProduct(A, epsV, epsU, xV)

            epsV_vals = epsV.evaluate(xV);
            epsU_vals = epsU.evaluate(xV);

            v11 = squeeze(epsV_vals(1,1,:,:));
            v22 = squeeze(epsV_vals(2,2,:,:));
            v12 = squeeze(epsV_vals(1,2,:,:));

            u11 = squeeze(epsU_vals(1,1,:,:));
            u22 = squeeze(epsU_vals(2,2,:,:));
            u12 = squeeze(epsU_vals(1,2,:,:));

            nPoints = numel(v11);
            vVoigt = [v11(:)'; v22(:)'; v12(:)'];
            uVoigt = [u11(:)'; u22(:)'; u12(:)'];

            temp = A * uVoigt;
            result = sum(vVoigt .* temp, 1);

            % Reshape
            result = reshape(result, size(v11));
            result = permute(result, [3, 1, 2]);
        end

        %% customDP
        function df = customDP(H, v, gradU)

            s.operation = @(xV) TutorialShells.evaluateCustomDP(H, v, gradU, xV);
            s.mesh = v.mesh;
            s.ndimf = 1;
            df = DomainFunction(s);
        end

        function result = evaluateCustomDP(H, v, gradU, xV)
            v_vals = v.evaluate(xV);       % 2 x nGauss x nElem
            gradU_vals = gradU.evaluate(xV); % 1 x 2 x nGauss x nElem (dimensión extra!)

            
            v1 = squeeze(v_vals(1,:,:));     % nGauss x nElem
            v2 = squeeze(v_vals(2,:,:));     % nGauss x nElem

            gu1 = squeeze(gradU_vals(1,1,:,:)); % nGauss x nElem
            gu2 = squeeze(gradU_vals(1,2,:,:)); % nGauss x nElem

            nPoints = numel(v1);
            vVec = [v1(:)'; v2(:)'];           % 2 x nPoints
            guVec = [gu1(:)'; gu2(:)'];        % 2 x nPoints

            temp = H * guVec;  % 2 x nPoints
            result = sum(vVec .* temp, 1);  % 1 x nPoints

            % Reshape
            result = reshape(result, size(v1));
            result = permute(result, [3, 1, 2]);  % 1 x nGauss x nElem
        end

        %% customDPgrad
        function df = customDPgrad(H, gradV, gradU)

            s.operation = @(xV) TutorialShells.evaluateCustomDPgrad(H, gradV, gradU, xV);
            s.mesh = gradV.mesh;
            s.ndimf = 1;
            df = DomainFunction(s);
        end

        function result = evaluateCustomDPgrad(H, gradV, gradU, xV)
            gradV_vals = gradV.evaluate(xV);
            gradU_vals = gradU.evaluate(xV);

            % Formato (1, 2, nGauss, nElem)
            gv1 = squeeze(gradV_vals(1,1,:,:));
            gv2 = squeeze(gradV_vals(1,2,:,:));
            gu1 = squeeze(gradU_vals(1,1,:,:));
            gu2 = squeeze(gradU_vals(1,2,:,:));

            nPoints = numel(gv1);
            gvVec = [gv1(:)'; gv2(:)'];
            guVec = [gu1(:)'; gu2(:)'];

            temp = H * guVec;
            result = sum(gvVec .* temp, 1);

            result = reshape(result, size(gv1));
            result = permute(result, [3, 1, 2]);
        end

        %% customDPvector
        function df = customDPvector(H, v, u)

            s.operation = @(xV) TutorialShells.evaluateCustomDPvector(H, v, u, xV);
            s.mesh = v.mesh;
            s.ndimf = 1;
            df = DomainFunction(s);
        end

        function result = evaluateCustomDPvector(H, v, u, xV)
            v_vals = v.evaluate(xV);  % 2 x nGauss x nElem
            u_vals = u.evaluate(xV);  % 2 x nGauss x nElem

            v1 = squeeze(v_vals(1,:,:));
            v2 = squeeze(v_vals(2,:,:));
            u1 = squeeze(u_vals(1,:,:));
            u2 = squeeze(u_vals(2,:,:));

            nPoints = numel(v1);
            vVec = [v1(:)'; v2(:)'];
            uVec = [u1(:)'; u2(:)'];

            temp = H * uVec;
            result = sum(vVec .* temp, 1);

            result = reshape(result, size(v1));
            result = permute(result, [3, 1, 2]);
        end

    end

end