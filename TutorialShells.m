classdef TutorialShells < handle


    properties (Access=private)
    % properties (Access = {?TutorialShells, ?auxiliaryFunctions})  %(Access = private)
        mesh
        young
        area
        shear
        inertia
        poisson, rotation 
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
        materialProperties, stressCase
    end

    methods (Access = public)
        %% TutorialShells
        function obj = TutorialShells()
            close all; 
            
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
            [strainFun, stressFun] = obj.createStrainStressFunctions(z, epsilonU_nodal, epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal);

            % Plane Stresses: obtain Strain and Stresses Functions on a
            % certain height for epsilon_zz = sigma_zz = 0
            [planeStrain, planeStress] = obj.planeStrainStress(h, epsilonU_nodal, epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal);
            % 
            % % Obtain Internal Forces 
            % [Nfun, Mfun, Qfun] = obj.internalForces(epsilonU_nodal,epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal);

            %% 
            [maxW, significantValues] = obj.findMaxDisplacements();
            


           

            %% PLOT AND PRINT   

            plotMatlab = false;
            printParaview = true;

            % KAPPA ==========================================
            kappa = 2;

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

                % % Plane Strains
                % obj.customPlot(planeStrain{kappa}, {'\epsilon_{xx}', '\epsilon_{yy}', '\epsilon_{xy}', '\epsilon_{xz}', '\epsilon_{yz}'});
                % 
                % % Plane Stresses 
                % obj.customPlot(planeStress{kappa}, {'\sigma_{xx}', '\sigma_{yy}', '\sigma_{xy}','\sigma_{xz}','\sigma_{yz}'});


                % % Resulting forcees 
                % obj.customPlot(Nfun, {'N_{xx}', 'N_{yy}', 'N_{xy}'});
                % obj.customPlot(Mfun, {'M_{xx}', 'M_{yy}', 'M_{xy}'});
                % obj.customPlot(Qfun, {'Q_{xz}', 'Q_{yz}'});
            end


            if printParaview == true
                obj.wFun.print('wfun print','Paraview')
                % obj.uFun.print('ufun print','Paraview')
                % obj.thetaFun.print('thetafun print','Paraview')
                % strainFun{kappa}.print('strain', 'Paraview'); % Kappa defined on plots
                % stressFun{kappa}.print('stress', 'Paraview');
                
                % stressFun{1}.print('BOTTOM Layer 1','Paraview')
                % stressFun{2}.print('TOP Layer 1','Paraview')
                % stressFun{3}.print('BOTTOM Layer 2','Paraview')
                % stressFun{4}.print('TOP Layer 2','Paraview')

                % Nfun.print('Nfun','Paraview');
                % Mfun.print('Mfun','Paraview');
                % Qfun.print('Qfun','Paraview');


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
            % -------------------------------------------------------------------------
            % MATERIAL PROPERTIES INPUT FORMAT
            % -------------------------------------------------------------------------
            %
            % This function allows defining material properties for:
            %   - ISOTROPIC materials
            %   - ORTHOTROPIC (orthotropic lamina) materials
            %
            % The material can be SINGLE layer or MULTI-layer laminate.
            %
            % -------------------------------------------------------------------------
            % -------------------------- ISOTROPIC CASE -------------------------------
            % -------------------------------------------------------------------------
            %
            % Set:
            %   obj.materialProperties = 'ISOTROPIC'
            %
            % Required inputs:
            %
            %   E  = [E1; E2; ...; En]          % Young modulus per layer
            %   nu = [nu1; nu2; ...; nun]       % Poisson ratio per layer
            %   h  = [h1; h2; ...; hn]          % Thickness of each layer
            %
            % Shear modulus is computed automatically as:
            %
            %   G = E ./ (2*(1+nu))
            %
            % Notes:
            % - Each row corresponds to one layer.
            % - For isotropic materials, properties are identical in all directions.
            % - Only one Young modulus and one Poisson ratio per layer are required.
            %
            %
            % -------------------------------------------------------------------------
            % ------------------------ ORTHOTROPIC CASE -------------------------------
            % -------------------------------------------------------------------------
            %
            % Set:
            %   obj.materialProperties = 'ORTHOTROPIC'
            %
            % Required inputs per layer:
            %
            %   E  = [E1 E2 E3]                % Young moduli in material directions
            %                                    E1 = longitudinal
            %                                    E2 = transverse
            %                                    E3 = through-thickness
            %
            %   nu = [nu12 nu13 nu23]          % Poisson ratios
            %
            %   G  = [G12 G13 G23]             % Shear moduli
            %
            %   h  = [h1; h2; ...; hn]         % Thickness of each layer
            %
            %   Rotation = [theta1; theta2; ...; thetan]  % Ply orientation in degrees
            %
            % Notes:
            % - Directions (1,2,3) are material principal axes.
            % - Rotation defines fibre orientation with respect to global axes.
            % - If MULTI-layer laminate is used:
            %       length(h) must equal number of plies
            %       length(Rotation) must match number of plies
            %
            % -------------------------------------------------------------------------
            % OUTPUT
            % -------------------------------------------------------------------------
            %
            % The function automatically computes:
            %
            %   A_Matrix  -> Membrane stiffness matrix
            %   B_Matrix  -> Membrane-bending coupling matrix
            %   D_Matrix  -> Bending stiffness matrix
            %   H_Matrix  -> Shear stiffness matrix
            %
            % using Classical Laminate Theory integration through thickness.
            %
            % -------------------------------------------------------------------------

            % --------------- DATABASE INPUT -------------------
            % Available isotropic materials:
            %   'Aluminum' - Aluminum alloy (E=10.6 msi, nu=0.33)
            %   'Copper'   - Pure copper (E=18.0 msi, nu=0.33)
            %   'Steel'    - Structural steel (E=30.0 msi, nu=0.29)
            % Available orthotropic composite materials:
            %   'AS'   - Graphite-Epoxy AS/3501 (High-strength carbon fiber)
            %   'EpT'  - Graphite-Epoxy T300/934 (Standard-modulus carbon fiber)
            %   'Ep1'  - Glass-Epoxy type 1 (E-glass, high strength)
            %   'Ep2'  - Glass-Epoxy type 2 (E-glass, lower modulus)
            %   'BrEp' - Boron-Epoxy (High stiffness, aerospace grade)
            %
            % --------------- INPUT -----------------------------
            % materialName = {'Material1'; 'Material2'; 'Material3'}

            materialName = {'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'};
            max_thickness = 0.5;
            h = max_thickness/length(materialName)*ones(length(materialName),1);

            % Rotation = [0, 0, 45, 90, -45, 90, 45, -45, 0, 0].';                          % degrees
            % Rotation = [0, 0, 45, 90, -45, 45, 90, -45, 0, 0].';                          % degrees
            % Rotation = [0, 0, 45, -45, 90, 90, -45, 45, 0, 0].';                          % degrees
            Rotation = [0, 0, 45, -45, 90, 90, -45, 45, 0, 0].';                          % degrees




            [E, nu, G, materialType] = obj.getMaterialProperties(materialName);
            obj.materialProperties = materialType;


            % --------------- MANUAL INPUT ---------------------
            % obj.materialProperties = 'ISOTROPIC'; % ISOTROPIC / ORTHOTROPIC
            %
            % % Isotropic 
            % E = [10; 3; 10];
            % nu = [0.3; 0.25; 0.3]; 
            % h = [0.25; 0.5; 0.25];
            % G = E ./ 2./ (1+nu);
            % 
            % % AS
            % E = [20 1.3 1.3]*1e6*6894.76;
            % nu = [0.3 0.3 0.49];                     % nu_12, nu_13, nu_23
            % h = 1;
            % G = [1.03 1.03 0.9]*1e6*6894.76;         % G_12, G_13, G_23
            % Rotation = [0];                          % degrees 

            

            % Set-up 
            obj.young = ConstantFunction.create(E,obj.mesh);
            obj.area = ConstantFunction.create(h,obj.mesh);  % CONSIDERED AS THICKNESS
            obj.shear = ConstantFunction.create(G,obj.mesh);
            obj.inertia = ConstantFunction.create(1,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);
            switch obj.materialProperties
                case 'ORTHOTROPIC'
                    obj.rotation = ConstantFunction.create(deg2rad(Rotation),obj);
            end

            % Layer properties (SCALARS)
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

            obj.stressCase = 'NORMAL';

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
            obj.bcCase = 1;     % 1 --> Modified bc // 2 --> Original bc (change q) 
            
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

                            pointloadFun = [];
                            for i = 1:numel(sPL)
                                pl = TractionLoad2(obj.mesh, sPL{i}, 'DIRAC');
                                pointloadFun = [pointloadFun, pl];

                            end
                            s.pointloadFun = pointloadFun;

                        case 2

                            % % Load on Right edge // CHANGE: q = 100
                            % sPL{1}.domain    = @(coor) isRight(coor);
                            % sPL{1}.direction = 2;
                            % sPL{1}.value     = 1;

                            s.pointloadFun = [];
                    end

                    

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
                    q = ConstantFunction.create(100,obj.mesh);
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

            matrixCase = 2; 
            % 1 = Scalar value 
            % 2 = Matrix value 

            switch matrixCase 
                case 1 
                    A = obj.A_coupling;
                    f = @(u,v) A.*DDP(SymGrad(v),SymGrad(u));
                case 2
                    A = obj.A_Matrix;
                    f = @(u,v) TutorialShells.customDDP(A, SymGrad(v), SymGrad(u));
            end

            Ku = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Ku = reduceMatrix(obj,Ku,obj.bcU,obj.bcU);

            

            switch matrixCase
                case 1
                    D = obj.D_bending;
                    f = @(u,v) D.*DDP(SymGrad(v),SymGrad(u));
                case 2
                    D = obj.D_Matrix;
                    f = @(u,v) TutorialShells.customDDP(D, SymGrad(v), SymGrad(u));
            end

            Ktheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Ktheta = obj.reduceMatrix(Ktheta,obj.bcT,obj.bcT);

            switch matrixCase
                case 1
                    B = obj.B_coupling;
                    f = @ (u,v) B.*DDP(SymGrad(v),SymGrad(u));
                case 2 
                    B = obj.B_Matrix;
                    f = @ (u,v) TutorialShells.customDDP(B,SymGrad(v),SymGrad(u));
            end

            Zut = IntegrateLHS(f,obj.uFun,obj.thetaFun,obj.mesh,'Domain',2);
            Zut = obj.reduceMatrix(Zut,obj.bcU,obj.bcT);

            switch matrixCase
                case 1
                    H = obj.H_shear;
                    f = @(u,v) H.*DP(v,u);
                case 2
                    H = obj.H_Matrix;
                    f = @(u,v) TutorialShells.customDPvector(H, v, u);
            end

            Mtheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Mtheta = obj.reduceMatrix(Mtheta,obj.bcT,obj.bcT);

            switch matrixCase
                case 1
                    f = @(u,v) H.*DP(v,Grad(u));
                case 2
                    f = @(u,v) TutorialShells.customDP(H, v, Grad(u));
            end

            Nthetaw = IntegrateLHS(f,obj.thetaFun,obj.wFun,obj.mesh,'Domain',2);            
            Nthetaw = obj.reduceMatrix(Nthetaw,obj.bcT,obj.bcW);
            
            switch matrixCase
                case 1
                    f = @(u,v) H.*DP(Grad(v),Grad(u));
                case 2
                    f = @(u,v) TutorialShells.customDPgrad(H, Grad(v), Grad(u));
            end

            Kw = IntegrateLHS(f,obj.wFun,obj.wFun,obj.mesh,'Domain',2);  
            Kw = obj.reduceMatrix(Kw,obj.bcW,obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            beta = 1; 

            % Zut = zeros(nU,nTheta);
            Zuw = zeros(nU,nW);
            LHS = [Ku Zut Zuw; Zut' (Ktheta+beta*Mtheta) beta*Nthetaw; Zuw' beta*Nthetaw' beta*Kw];
        end

        %% createStrainStressFunctions
        function [strainFun, stressFun] = createStrainStressFunctions(obj, z, epsilonU_nodal, epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal)

            fprintf('\n===== COMPLETE STRESSES =====\n');
            
            obj.stressCase = 'NORMAL';

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
            switch obj.stressCase
                case 'NORMAL'
                    kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);
                    % when there is layer superposition, takes the value of the bottom layer

                    switch obj.materialProperties
                        case 'ISOTROPIC'
                            % kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);
                            % when there is layer superposition, takes the value of the bottom layer

                            if isInterface == false
                                nu = obj.poisson.constant(kLayer);
                                lambda = obj.young.constant(kLayer) / (1+nu) / (1-2*nu);
                               
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
                                    
                                    C(:,:,i) = lambda* [1-nu, nu, nu, 0, 0, 0;
                                        nu, 1-nu, nu, 0, 0, 0;
                                        nu, nu, 1-nu, 0, 0, 0;
                                        0, 0, 0, (1-2*nu), 0, 0;
                                        0, 0, 0, 0, (1-2*nu), 0;
                                        0, 0, 0, 0, 0, (1-2*nu)];

                                    kLayer = kLayer + 1;
                                end
                            end
                        case 'ORTHOTROPIC'
                            % kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);

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

                                theta = obj.rotation.constant(kLayer,1);

                                % T = [cos(theta)^2,         sin(theta)^2,          0,  0,           0,          -sin(2*theta);
                                %     sin(theta)^2,          cos(theta)^2,          0,  0,           0,           sin(2*theta);
                                %     0,                     0,                     1,  0,           0,           0;
                                %     sin(theta)*cos(theta), -sin(theta)*cos(theta), 0,  0,           0,           cos(theta)^2 - sin(theta)^2
                                %     0,                     0,                     0, -sin(theta),  cos(theta),  0;
                                %     0,                     0,                     0,  cos(theta),  sin(theta),  0];
                                T = [cos(theta)^2,         sin(theta)^2,          0, -sin(2*theta),  0,           0;
                                        sin(theta)^2,          cos(theta)^2,          0, sin(2*theta),  0,           0;
                                        0,                     0,                     1,  0,           0,           0;
                                        sin(theta)*cos(theta), -sin(theta)*cos(theta), 0, cos(theta)^2 - sin(theta)^2,  0,           0;
                                        0,                     0,                     0, 0,  cos(theta),  -sin(theta);
                                        0,                     0,                     0,  0,  sin(theta),  cos(theta)];

                                C = T* [ C11 C12 C13 0   0   0;
                                    C12 C22 C23 0   0   0;
                                    C13 C23 C33 0   0   0;
                                    0   0   0   2*C44 0   0;
                                    0   0   0   0   2*C55 0;
                                    0   0   0   0   0   2*C66 ]*T.';
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

                                    theta = obj.rotation.constant(kLayer,1);

                                    T = [cos(theta)^2,         sin(theta)^2,          0, -sin(2*theta),  0,           0;
                                        sin(theta)^2,          cos(theta)^2,          0, sin(2*theta),  0,           0;
                                        0,                     0,                     1,  0,           0,           0;
                                        sin(theta)*cos(theta), -sin(theta)*cos(theta), 0, cos(theta)^2 - sin(theta)^2,  0,           0;
                                        0,                     0,                     0, 0,  cos(theta),  -sin(theta);
                                        0,                     0,                     0,  0,  sin(theta),  cos(theta)];

                                    C(:,:,i) = T*[ C11 C12 C13 0   0   0;
                                        C12 C22 C23 0   0   0;
                                        C13 C23 C33 0   0   0;
                                        0   0   0   2*C44 0   0;
                                        0   0   0   0   2*C55 0;
                                        0   0   0   0   0   2*C66 ]*T.';

                                    kLayer = kLayer+1;
                                end
                            end
                    end
                case 'PLANE'
                    kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);
                    % when there is layer superposition, takes the value of the bottom layer
                    switch obj.materialProperties
                        case 'ISOTROPIC'
                            % kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);
                            % when there is layer superposition, takes the value of the bottom layer

                            if isInterface == false
                                nu = obj.poisson.constant(kLayer);
                                lambda = obj.young.constant(kLayer) / (1+nu) / (1-2*nu);
                               
                                C = lambda* [1-nu, nu, 0, 0, 0;
                                    nu, 1-nu, 0, 0, 0;
                                    0, 0, (1-2*nu), 0, 0;
                                    0, 0, 0, (1-2*nu), 0;
                                    0, 0, 0, 0, (1-2*nu)];
                            else
                                for i = 1:2
                                    nu = obj.poisson.constant(kLayer);
                                    lambda = obj.young.constant(kLayer) / (1+nu) / (1-2*nu);
                                    
                                    C(:,:,i) = lambda* [1-nu, nu, 0, 0, 0;
                                        nu, 1-nu, 0, 0, 0;
                                        0, 0, (1-2*nu), 0, 0;
                                        0, 0, 0, (1-2*nu), 0;
                                        0, 0, 0, 0, (1-2*nu)];

                                    kLayer = kLayer + 1;
                                end
                            end
                        case 'ORTHOTROPIC'
                            % kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);

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

                                theta = obj.rotation.constant(kLayer,1);

                                T = [cos(theta)^2,         sin(theta)^2,          0,           0,          -sin(2*theta);
                                    sin(theta)^2,          cos(theta)^2,          0,           0,           sin(2*theta);
                                    0,                     0,                     cos(theta),  sin(theta),  0;
                                    0,                     0,                     -sin(theta),  cos(theta),  0;
                                    sin(theta)*cos(theta), -sin(theta)*cos(theta), 0,           0,           cos(theta)^2 - sin(theta)^2];

                                C = T* [ C11 C12  0   0   0;
                                    C12 C22  0   0   0;
                                    0   0   2*C44 0   0;
                                    0   0   0   2*C55 0;
                                    0   0   0   0   2*C66 ]*T.';
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
                                    C22 = (1 - nu13*nu31)/(E1*E3*delta);
                                    C44 = G12;
                                    C55 = G13;
                                    C66 = G23;

                                    theta = obj.rotation.constant(kLayer,1);

                                    T = [cos(theta)^2,         sin(theta)^2,          0,           0,          -sin(2*theta);
                                        sin(theta)^2,          cos(theta)^2,          0,           0,           sin(2*theta);
                                        0,                     0,                     cos(theta),  sin(theta),  0;
                                        0,                     0,                     -sin(theta),  cos(theta),  0;
                                        sin(theta)*cos(theta), -sin(theta)*cos(theta),0,           0,           cos(theta)^2 - sin(theta)^2];

                                    C(:,:,i) = T*[ C11 C12 0   0   0;
                                        C12 C22 0   0   0;
                                        0   0   2*C44 0   0;
                                        0   0   0   2*C55 0;
                                        0   0   0   0   2*C66 ]*T.';

                                    kLayer = kLayer+1;
                                end
                            end
                    end

            end
        end

        %% planeStrainStress

        function [strainFun, stressFun] = planeStrainStress(obj,z, epsilonU_nodal, epsilonTheta_nodal, dw_dx_nodal, dw_dy_nodal)

            fprintf('\n===== PLANE STRESSES =====\n');

            obj.stressCase = 'PLANE';

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
                strain_nodal = zeros(nNodes, 5);
                strain_nodal(:, 1) = epsilonU_nodal(1, :) + z{i} * epsilonTheta_nodal(1, :);  % exx
                strain_nodal(:, 2) = epsilonU_nodal(2, :) + z{i} * epsilonTheta_nodal(2, :);  % eyy
                strain_nodal(:, 3) = epsilonU_nodal(3, :) + z{i} * epsilonTheta_nodal(3, :);  % exy
                strain_nodal(:, 4) = epsilon_xz_nodal;
                strain_nodal(:, 5) = epsilon_yz_nodal;

                % Create LagrangianFunction for strains [exx, eyy, ezz, exy, exz, eyz]
                strainFun{i} = LagrangianFunction.create(obj.mesh, 5, 'P1');
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
                    stressTemp{1,i} = LagrangianFunction.create(obj.mesh, 5, 'P1');
                    stressTemp{1,i}.setFValues(stress_nodal);
                else 
                    for j = 1:2
                        stress_nodal = (C(:,:,j) * strain_nodal.').';  % (nNodes x 6)
                        stressTemp{j,i} = LagrangianFunction.create(obj.mesh, 5, 'P1');
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

        %% getMaterialProperties
        function [E, nu, G, type] = getMaterialProperties(obj,materialName)
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

            db.Copper.type = 'ISOTROPIC';
            db.Copper.E  = 18.0 * msi_to_Pa;
            db.Copper.nu = 0.33;
            db.Copper.G  = 6.39 * msi_to_Pa;

            db.Steel.type = 'ISOTROPIC';
            db.Steel.E  = 30.0 * msi_to_Pa;
            db.Steel.nu = 0.29;
            db.Steel.G  = 11.24 * msi_to_Pa;

            % ORTHOTROPIC
            db.AS.type = 'ORTHOTROPIC';
            db.AS.E  = [20.0, 1.3, 1.3] * msi_to_Pa;
            db.AS.nu = [0.30, 0.30, 0.49];
            db.AS.G  = [1.03, 1.03, 0.90] * msi_to_Pa;

            db.EpT.type = 'ORTHOTROPIC';
            db.EpT.E  = [19.0, 1.5, 1.5] * msi_to_Pa;
            db.EpT.nu = [0.22, 0.22, 0.49];
            db.EpT.G  = [1.00, 0.90, 0.90] * msi_to_Pa;

            db.Ep1.type = 'ORTHOTROPIC';
            db.Ep1.E  = [7.8, 2.6, 2.6] * msi_to_Pa;
            db.Ep1.nu = [0.25, 0.25, 0.34];
            db.Ep1.G  = [1.30, 1.30, 0.50] * msi_to_Pa;

            db.Ep2.type = 'ORTHOTROPIC';
            db.Ep2.E  = [5.6, 1.2, 1.3] * msi_to_Pa;
            db.Ep2.nu = [0.26, 0.26, 0.34];
            db.Ep2.G  = [0.60, 0.60, 0.50] * msi_to_Pa;

            db.BrEp.type = 'ORTHOTROPIC';
            db.BrEp.E  = [30.0, 3.0, 3.0] * msi_to_Pa;
            db.BrEp.nu = [0.30, 0.25, 0.25];
            db.BrEp.G  = [1.00, 1.00, 0.60] * msi_to_Pa;

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

                for i = 1:nMaterials
                    E(i)  = db.(fieldNames{i}).E;
                    nu(i) = db.(fieldNames{i}).nu;
                    G(i)  = db.(fieldNames{i}).G;
                end
            else
                E  = zeros(nMaterials, 3);
                nu = zeros(nMaterials, 3);
                G  = zeros(nMaterials, 3);

                for i = 1:nMaterials
                    mat = db.(fieldNames{i});
                    if strcmp(mat.type, 'ISOTROPIC')
                        E(i,:)  = [mat.E, mat.E, mat.E];
                        nu(i,:) = [mat.nu, mat.nu, mat.nu];
                        G(i,:)  = [mat.G, mat.G, mat.G];
                    else
                        E(i,:)  = mat.E;
                        nu(i,:) = mat.nu;
                        G(i,:)  = mat.G;
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
                end
            end
        end

        %% findMaxDisplacements
        function [maxW, significantValues] = findMaxDisplacements(obj, nValues)
            % findMaxDisplacements - Find maximum and other significant large w displacements
            %
            % Usage:
            %   [maxW, significantValues] = findMaxDisplacements(obj)
            %   [maxW, significantValues] = findMaxDisplacements(obj, nValues)
            %
            % Inputs:
            %   obj     - TutorialShells object with wFun property
            %   nValues - (optional) Number of significant values to show (default: 10)
            %
            % Outputs:
            %   maxW              - Absolute maximum displacement value
            %   significantValues - Struct array with largest displacement values:
            %                       .value   - Displacement value
            %                       .abs     - Absolute value
            %                       .indices - Node indices where it occurs
            %                       .coords  - [x, y] coordinates
            %                       .count   - Number of nodes with this value

            if nargin < 2
                nValues = 10;  % Show top 10 by default
            end

            % Get w displacements
            w_values = obj.wFun.fValues;
            coords = obj.mesh.coord;
            nNodes = length(w_values);

            % Find absolute maximum
            [maxW_abs, idx_max_abs] = max(abs(w_values));
            maxW = w_values(idx_max_abs);

            fprintf('\n========== MAXIMUM DISPLACEMENT ==========\n');
            fprintf('Maximum |w| = %.6e\n', maxW_abs);
            fprintf('Value w    = %.6e\n', maxW);
            fprintf('Node index = %d\n', idx_max_abs);
            fprintf('Location   = (%.4f, %.4f)\n', coords(idx_max_abs, 1), coords(idx_max_abs, 2));

            % Sort all displacements by absolute value
            [w_abs_sorted, sort_idx] = sort(abs(w_values), 'descend');
            w_sorted = w_values(sort_idx);

            % Group similar values together (tolerance for floating point)
            tolerance = maxW_abs * 1e-10;  % Very small tolerance for grouping

            significantValues = struct('value', {}, 'abs', {}, 'indices', {}, 'coords', {}, 'count', {});

            fprintf('\n========== SIGNIFICANT LARGE DISPLACEMENTS ==========\n');
            fprintf('Rank |      Value       |   |Value|   |  %% of Max | Nodes | Locations\n');
            fprintf('-----+------------------+-------------+----------+-------+------------------------\n');

            idx = 1;
            rank = 1;
            processed = false(nNodes, 1);

            while rank <= nValues && idx <= nNodes
                if processed(sort_idx(idx))
                    idx = idx + 1;
                    continue;
                end

                current_val = w_sorted(idx);
                current_abs = w_abs_sorted(idx);

                % Find all nodes with similar value (within tolerance)
                similar_mask = abs(w_values - current_val) < tolerance;
                similar_idx = find(similar_mask);

                % Mark as processed
                processed(similar_idx) = true;

                % Store in output
                significantValues(rank).value = current_val;
                significantValues(rank).abs = current_abs;
                significantValues(rank).indices = similar_idx;
                significantValues(rank).coords = coords(similar_idx, :);
                significantValues(rank).count = length(similar_idx);

                % Print
                percent_of_max = (current_abs / maxW_abs) * 100;

                fprintf('%3d  | %+15.6e | %12.6e | %7.2f%% | %5d | ', ...
                    rank, current_val, current_abs, percent_of_max, length(similar_idx));

                % Show first few node coordinates
                if length(similar_idx) <= 3
                    for i = 1:length(similar_idx)
                        fprintf('(%4.2f,%4.2f) ', coords(similar_idx(i), 1), coords(similar_idx(i), 2));
                    end
                else
                    fprintf('(%4.2f,%4.2f) (%4.2f,%4.2f) ... +%d more', ...
                        coords(similar_idx(1), 1), coords(similar_idx(1), 2), ...
                        coords(similar_idx(2), 1), coords(similar_idx(2), 2), ...
                        length(similar_idx) - 2);
                end
                fprintf('\n');

                rank = rank + 1;
                idx = idx + 1;
            end

            fprintf('================================================\n');

            % Statistical summary
            fprintf('\n========== STATISTICAL SUMMARY ==========\n');
            fprintf('Total nodes:           %d\n', nNodes);
            fprintf('Max |w|:               %.6e\n', maxW_abs);
            fprintf('Min |w|:               %.6e\n', min(abs(w_values)));
            fprintf('Mean |w|:              %.6e\n', mean(abs(w_values)));
            fprintf('Std |w|:               %.6e\n', std(abs(w_values)));
            fprintf('Nodes > 50%% max:       %d (%.1f%%)\n', sum(abs(w_values) > 0.5*maxW_abs), ...
                100*sum(abs(w_values) > 0.5*maxW_abs)/nNodes);
            fprintf('Nodes > 10%% max:       %d (%.1f%%)\n', sum(abs(w_values) > 0.1*maxW_abs), ...
                100*sum(abs(w_values) > 0.1*maxW_abs)/nNodes);
            fprintf('=========================================\n\n');

            % % Visualization
            % figure('Name', 'Maximum Displacement Analysis', 'Position', [100, 100, 1200, 500]);
            % 
            % % Plot 1: Full displacement field with top locations marked
            % subplot(1,2,1);
            % trisurf(obj.mesh.connec, coords(:,1), coords(:,2), w_values);
            % view(0, 90);
            % shading interp;
            % colorbar;
            % hold on;
            % 
            % % Mark top significant values with different colors
            % colors = lines(min(5, length(significantValues)));
            % for i = 1:min(5, length(significantValues))
            %     idx_sig = significantValues(i).indices;
            %     plot3(coords(idx_sig, 1), coords(idx_sig, 2), ...
            %         w_values(idx_sig), 'o', 'Color', colors(i,:), ...
            %         'MarkerSize', 10, 'LineWidth', 2, ...
            %         'DisplayName', sprintf('Rank %d: %.2e', i, significantValues(i).value));
            % end
            % 
            % legend('show', 'Location', 'best');
            % title('w Displacement Field with Significant Values');
            % xlabel('x'); ylabel('y');
            % grid on;
            % 
            % % Plot 2: Absolute values
            % subplot(1,2,2);
            % trisurf(obj.mesh.connec, coords(:,1), coords(:,2), abs(w_values));
            % view(0, 90);
            % shading interp;
            % colorbar;
            % hold on;
            % 
            % % Mark maximum
            % plot3(coords(idx_max_abs, 1), coords(idx_max_abs, 2), maxW_abs, ...
            %     'r*', 'MarkerSize', 20, 'LineWidth', 3, 'DisplayName', 'Maximum');
            % 
            % legend('show');
            % title('|w| Displacement Field');
            % xlabel('x'); ylabel('y');
            % grid on;
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