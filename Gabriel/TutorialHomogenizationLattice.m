classdef TutorialHomogenizationLattice < handle

    properties (Access = public)
        paramHole
        Chomog
        volFrac
        f,df,ddf
    end

    properties (Access = private)
        E
        nu
        meshType
        meshN
        holeType
        nSteps
        % damageType
        pnorm
        monitoring
        Mmass
        fixedHoleSize
        currentB
        latticeVectors
    end

    properties (Access = private)
        baseMesh
        masterSlave
        test
        maxParam
        
    end

    methods (Access = public)
        
        function obj = TutorialHomogenizationLattice()
            obj.init();
            % obj.defineMesh();
            obj.computeHoleParams();
            obj.compute();
            obj.fitting();
            obj.plot();
            obj.printTensorAtVolume(0.5);
        end 
        function compareFunctionsAofB(obj)
            % Compara 5 funções a(b) e prova que nenhuma simétrica
            % satisfaz a(0)=1 E a<1 E a>1 simultaneamente
            
            b_vals = linspace(-0.6, 0.6, obj.nSteps);
            
            % Definir as funções
            funcs = {
                @(b) exp(b.^2) - 0.5,           'exp(b²)-0.5  [a(0)=0.5 ✗]';
                @(b) exp(b.^2 * log(1.5)),       'exp(b²·ln1.5) [a(0)=1 ✓ mas a≥1]';
                @(b) 1 + 0.6*b + 0.2*b.^2,      '1+0.6b+0.2b² [assimétrica]';
                @(b) exp(b.^2),                  'exp(b²)  [a(0)=1 ✓ mas a≥1]';
                @(b) 1./(1 + 0.5*b.^2),         '1/(1+0.5b²) [a(0)=1 ✓ mas a≤1]';
            };
            
            nFuncs = size(funcs, 1);
            colors = lines(nFuncs);
            
            % ===== Plot 1: a(b) para cada função =====
            figure('Position', [50, 50, 1800, 1000], 'Color', 'white');
            
            subplot(2, nFuncs+1, 1);
            hold on;
            for k = 1:nFuncs
                a_vals = funcs{k,1}(b_vals);
                plot(b_vals, a_vals, '-', 'Color', colors(k,:), 'LineWidth', 2, ...
                     'DisplayName', funcs{k,2});
            end
            yline(1, 'k--', 'a=1', 'LineWidth', 1.5);
            xlabel('b'); ylabel('a(b)'); title('Funções a(b)'); grid on;
            legend('Location', 'best', 'FontSize', 7);
            
            % ===== Para cada função: calcular C^h e plotar =====
            for k = 1:nFuncs
                fprintf('\n=== Funcion %d: %s ===\n', k, funcs{k,2});
                
                C11_vals = zeros(1, obj.nSteps);
                C22_vals = zeros(1, obj.nSteps);
                C33_vals = zeros(1, obj.nSteps);
                meshes   = cell(1, obj.nSteps);
                
                for i = 1:obj.nSteps
                    b_val = b_vals(i);
                    a_val = funcs{k,1}(b_val);
                    
                    % Proteger contra a muito pequeno
                    if a_val < 0.1
                        C11_vals(i) = NaN;
                        C22_vals(i) = NaN;
                        C33_vals(i) = NaN;
                        continue;
                    end
                    
                    d_val = (1 + b_val^2) / a_val;
                    obj.latticeVectors = [a_val, b_val; b_val, d_val];
                    obj.defineMesh();
                    
                    C = obj.computeHomogenization(b_val);
                    C11_vals(i) = C(1,1,1,1);
                    C22_vals(i) = C(2,2,2,2);
                    C33_vals(i) = C(1,2,1,2);
                    
                    % Guardar malhas para 3 pontos representativos
                    if i == 1 || i == round(obj.nSteps/2) || i == obj.nSteps
                        meshes{i} = obj.baseMesh;
                    end
                    
                    if mod(i, 10) == 0
                        fprintf('  %d/%d b=%.3f a=%.3f\n', i, obj.nSteps, b_val, a_val);
                    end
                end
                
                % Plot C11, C22, C33
                subplot(2, nFuncs+1, k+1);
                hold on;
                plot(b_vals, C11_vals, 'b-', 'LineWidth', 1.5, 'DisplayName', 'C_{11}');
                plot(b_vals, C22_vals, 'r-', 'LineWidth', 1.5, 'DisplayName', 'C_{22}');
                plot(b_vals, C33_vals, 'g-', 'LineWidth', 1.5, 'DisplayName', 'C_{33}');
                xlabel('b'); title(funcs{k,2}, 'FontSize', 8, 'Interpreter', 'none');
                legend('Location', 'best', 'FontSize', 7); grid on;
                
                % Plot 3 micros representativas
                idx_pts = [1, round(obj.nSteps/2), obj.nSteps];
                b_pts   = b_vals(idx_pts);
                
                for m = 1:3
                    subplot(2, nFuncs+1, (nFuncs+1) + k);
                    if m == 1; hold on; axis off; axis equal; end
                    
                    b_local = b_pts(m);
                    a_local = funcs{k,1}(b_local);
                    if a_local < 0.1; continue; end
                    d_local = (1 + b_local^2) / a_local;
                    
                    % Posição horizontal das 3 micros
                    offset_x = (m-1) * 2.5;
                    scale    = 0.8;
                    v1 = scale * [a_local, b_local];
                    v2 = scale * [b_local, d_local];
                    
                   centro = [offset_x + 0.5*(v1(1)+v2(1)), 0.5*(v1(2)+v2(2))];
    
                    % Vértices do paralelogramo (a partir do centro)
                    px = centro(1) + [0, v1(1), v1(1)+v2(1), v2(1), 0] - 0.5*(v1(1)+v2(1));
                    py = centro(2) + [0, v1(2), v1(2)+v2(2), v2(2), 0] - 0.5*(v1(2)+v2(2));
                    fill(px, py, colors(k,:), 'FaceAlpha', 0.4, ...
                         'EdgeColor', colors(k,:), 'LineWidth', 1.5);
                    
                    % Furo
                    hf = 0.3;
                    local_pts = [0.5-hf/2, 0.5-hf/2;
                                 0.5+hf/2, 0.5-hf/2;
                                 0.5+hf/2, 0.5+hf/2;
                                 0.5-hf/2, 0.5+hf/2];
                    hx = centro(1) + (local_pts(:,1)-0.5)*v1(1) + (local_pts(:,2)-0.5)*v2(1);
                    hy = centro(2) + (local_pts(:,1)-0.5)*v1(2) + (local_pts(:,2)-0.5)*v2(2);
                    fill(hx, hy, 'w', 'EdgeColor', 'none');
                    
                    text(offset_x + 0.4, -0.25, sprintf('b=%.2f\na=%.2f', b_local, a_local), ...
                         'FontSize', 7, 'HorizontalAlignment', 'center');
                end
                title(sprintf('Micros: %s', funcs{k,2}), 'FontSize', 7, 'Interpreter', 'none');
            end
            
            sgtitle({'Comparar a(b)', ...
                     'a(0)=1 E a<1 E a>1'}, ...
                    'FontSize', 11, 'FontWeight', 'bold');
        end
        
        function plotMultipleHoleSizes(obj, holeSizes)
            % Plota C11, C22, C33, anisotropia e Zener ratio
            
            if nargin < 2
                holeSizes = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
            end
            
            vf_values = 1 - holeSizes.^2;
            nHoles = length(holeSizes);
            colors = jet(nHoles);
            
            all_b = cell(nHoles, 1);
            all_C11 = cell(nHoles, 1);
            all_C22 = cell(nHoles, 1);
            all_C33 = cell(nHoles, 1);
            all_C12 = cell(nHoles, 1);
            
            originalHoleSize = obj.fixedHoleSize;
            
            for h = 1:nHoles
                fprintf('Processando vf = %.3f (l=%.2f)\n', vf_values(h), holeSizes(h));
                
                obj.fixedHoleSize = holeSizes(h);
                obj.computeHoleParams();
                obj.compute();
                obj.fitting();
                
                all_b{h} = obj.paramHole;
                all_C11{h} = squeeze(obj.Chomog(1,1,1,1,:));
                all_C22{h} = squeeze(obj.Chomog(2,2,2,2,:));
                all_C33{h} = squeeze(obj.Chomog(1,2,1,2,:));
                all_C12{h} = squeeze(obj.Chomog(1,1,2,2,:));
            end
            
            obj.fixedHoleSize = originalHoleSize;
            
            figure('Name', 'Comparação', 'Position', [50, 50, 1800, 800]);
            
            % Plot 1: C11
            subplot(2,3,1);
            hold on;
            for h = 1:nHoles
                plot(all_b{h}, all_C11{h}, '-', 'Color', colors(h,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('C_{11}'); title('C_{11} vs b'); grid on;
            
            % Plot 2: C22
            subplot(2,3,2);
            hold on;
            for h = 1:nHoles
                plot(all_b{h}, all_C22{h}, '-', 'Color', colors(h,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('C_{22}'); title('C_{22} vs b'); grid on;
            
            % Plot 3: C33
            subplot(2,3,3);
            hold on;
            for h = 1:nHoles
                plot(all_b{h}, all_C33{h}, '-', 'Color', colors(h,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('C_{33}'); title('C_{33} vs b'); grid on;
            
            % Plot 4: Anisotropia
            subplot(2,3,4);
            hold on;
            for h = 1:nHoles
                plot(all_b{h}, all_C11{h} ./ all_C22{h}, '-', 'Color', colors(h,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('C_{11}/C_{22}'); title('Anisotropia vs b'); grid on;
            yline(1, 'k--', 'Isotropia');
            
            % Plot 5: Zener
            subplot(2,3,5);
            hold on;
            for h = 1:nHoles
                Zener = 2 * all_C33{h} ./ (all_C11{h} - all_C12{h});
                plot(all_b{h}, Zener, '-', 'Color', colors(h,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('Zener Ratio'); title('Zener Ratio vs b'); grid on;
            yline(1, 'k--', 'Isotropia');
            
            % Plot 6: Legenda (simples)
            subplot(2,3,6);
            axis off;
            for h = 1:nHoles
                text(0.1, 0.95 - (h-1)*0.05, sprintf('  V_f = %.2f', vf_values(h)), ...
                     'Color', colors(h,:), 'FontSize', 10, 'VerticalAlignment', 'top');
            end
            title('Leyenda');
            
            sgtitle('Efecto del cambio del agujero/Vfrac');
        end
      
        

    end
    
    methods (Access = private)
        
        function init(obj)
            obj.E          = 1;
            obj.nu         = 0.3;
            obj.meshType   = 'Square';
            obj.meshN      = 150;

            obj.holeType   = 'Square';
            obj.pnorm      = 'Inf';
            % obj.damageType = 'Area';
            obj.nSteps     = 80;
            obj.monitoring = false;
            obj.fixedHoleSize = 0.7;
        end

        % function defineMesh(obj)
        %     switch obj.meshType
        %         case 'Square'
        %             s.c = [1.2,0.6];
        %             s.theta = [0,90];
        %             s.divUnit = obj.meshN;
        %             s.filename = '';
        %             MC = MeshCreator(s);
        %             MC.computeMeshNodes();
        % 
        %             % obj.lattice = MC.lattice;
        %         case 'Hexagon'
        %             s.c = [1.3,0.6,0.4];
        %             s.theta = [0,60,120];
        %             s.divUnit = obj.meshN;
        %             s.filename = '';
        %             MC = MeshCreator(s);
        %             MC.computeMeshNodes();
        %             % obj.lattice = MC.lattice;
        %     end
        %     s.coord  = MC.coord;
        %     s.connec = MC.connec;
        %     obj.baseMesh = Mesh.create(s);
        %     obj.masterSlave = MC.masterSlaveIndex;
        %     obj.test = LagrangianFunction.create(obj.baseMesh,1,'P1');
        %     obj.Mmass = IntegrateLHS(@(u,v) DP(v,u), obj.test, obj.test, obj.baseMesh, 'Domain');
        % end

        function defineMesh(obj)
            switch obj.meshType
                case 'Square'
                    s.latticeVectors = obj.latticeVectors;
                    s.divUnit = obj.meshN;
                    s.filename = '';
                    MC = MeshCreator(s);
                    MC.computeMeshNodes();
                    
                case 'Hexagon'
                    v1 = obj.latticeVectors(1, :);
                    v2 = obj.latticeVectors(2, :);
                    v3 = v2 - v1;
                    s.latticeVectors = [v1; v2; v3];
                    s.divUnit = obj.meshN;
                    s.filename = '';
                    MC = MeshCreator(s);
                    MC.computeMeshNodes();
            end
            
            s.coord  = MC.coord;
            s.connec = MC.connec;
            obj.baseMesh = Mesh.create(s);
            obj.masterSlave = MC.masterSlaveIndex;
            obj.test = LagrangianFunction.create(obj.baseMesh,1,'P1');
            obj.Mmass = IntegrateLHS(@(u,v) DP(v,u), obj.test, obj.test, obj.baseMesh, 'Domain');
        end

        % function computeHoleParams(obj)
        %     obj.maxParam = 0.979*ones(size(obj.nSteps));
        %     nParam = length(obj.maxParam);
        %     obj.paramHole = cell(1,nParam);
        %     for i=1:nParam
        %         obj.paramHole{i} = linspace(1e-9,obj.maxParam(i),obj.nSteps(i));
        %     end
        % end
        function computeHoleParams(obj)
            obj.maxParam = 0.6;
            obj.paramHole = linspace(0,obj.maxParam, obj.nSteps);
        end

        % function compute(obj)
        %     comb = table2array(combinations(obj.paramHole{:}));
        %     nComb = size(comb,1);
        %     mat = zeros(2,2,2,2,nComb);
        %     volF = zeros(1,nComb);
        %     for i=1:nComb
        %         hole = comb(i,:);
        %         if i==1
        %             hole = 1e-10*ones(size(hole));
        %         end
        %         mat(:,:,:,:,i) = obj.computeHomogenization(hole);
        %         volF(i)    = obj.computeVolumeFraction(hole);
        %     end
        %     obj.Chomog = obj.assembleResults(mat);
        %     obj.volFrac = obj.assembleResults(volF);
        % end
         function compute(obj)
            nComb = length(obj.paramHole);
            mat = zeros(2,2,2,2,nComb);
            volF = zeros(1,nComb);
            
           
            
            for i = 1:nComb
                b_val = obj.paramHole(i);
                obj.currentB = b_val;             
                
                a = exp(b_val^2);
                d = (1 + b_val^2) / a;
                
                v1 = [a,    b_val];
                v2 = [b_val, d   ];
                % v1 = [1, 0];
                % v2 = [b_val, 1];
                obj.latticeVectors = [v1; v2];                
                
                obj.defineMesh();                
               
                mat(:,:,:,:,i) = obj.computeHomogenization(b_val);                
                
                volF(i) = obj.computeVolumeFraction(b_val);
                
                if mod(i,5) == 0
                    fprintf('  %d/%d - b = %.4f - Volume = %.4f\n', i, nComb, b_val, volF(i));
                end
            end
            
            obj.Chomog = mat;
            obj.volFrac = volF;
        end

        
        function matHomog = computeHomogenization(obj, b_val)
            dens = obj.createDensityLevelSet(b_val); 
            mat  = obj.createDensityMaterial(dens);
            matHomog = obj.solveElasticMicroProblem(mat, dens);
        end

        function lsf = createDensityLevelSet(obj, b_val)
            
            ls = obj.computeLevelSet(obj.baseMesh, obj.fixedHoleSize);  
            sUm.backgroundMesh = obj.baseMesh;
            sUm.boundaryMesh   = obj.baseMesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);

            ls = CharacteristicFunction.create(uMesh);
            s.trial = obj.test;
            s.mesh = obj.baseMesh;
            f = FilterLump(s); 
            lsf = f.compute(ls, 2);
        end
        function ls = computeLevelSet(obj, mesh, l_fixo)
            gPar.type = obj.holeType;
            gPar.pnorm = obj.pnorm;   
            
            coord = mesh.coord;
            xmin = min(coord(:,1));
            xmax = max(coord(:,1));
            ymin = min(coord(:,2));
            ymax = max(coord(:,2));          
            center_x = (xmin + xmax)/2;
            center_y = (ymin + ymax)/2;
            gPar.xCoorCenter = center_x;
            gPar.yCoorCenter = center_y;
            
            
            if ~isempty(obj.latticeVectors)
                v1 = obj.latticeVectors(1, :);
                v2 = obj.latticeVectors(2, :);
                gPar.a1 = v1;
                gPar.a2 = v2;
                phi = atan2(v1(2), v1(1));
            else
                if size(coord,1) >= 4
                    v1 = coord(2,:) - coord(1,:);
                    v2 = coord(3,:) - coord(2,:);
                    gPar.a1 = v1;
                    gPar.a2 = v2;
                    phi = atan2(v1(2), v1(1));
                else
                    gPar.a1 = [1, 0];
                    gPar.a2 = [0, 1];
                    phi = 0;
                end
            end
            gPar.rotation = phi;
            
            
            switch obj.holeType
                case 'Circle'
                    gPar.radius = l_fixo/2;
                case 'Square'
                    gPar.length = l_fixo;
                case 'SmoothRectangle'
                    sx = l_fixo;
                    sy = l_fixo/2;           
                    gPar.xSide = sx;
                    gPar.ySide = sy;
                    gPar.pnorm = 16;
                case 'Ellipse'
                    gPar.type = "SmoothRectangle";
                    gPar.xSide  = l_fixo(1);
                    gPar.ySide  = l_fixo(2);
                    gPar.pnorm  = 2;  
                case 'SmoothHexagon'
                    gPar.radius = l_fixo;
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
                case 'ReinforcedHoneycomb'
                    gPar.theta  = 1-l_fixo;                          
                    gPar.eps    = 1;                        
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];            
                    gPar.radius = l_fixo;    
                    gPar.rotation = phi;
            end  
            
            g = GeometricalFunction(gPar);
            phiFun = g.computeLevelSetFunction(mesh);
            lsCircle = phiFun.fValues;
            ls = -lsCircle;            
        end
        
        % function ls = computeLevelSet(obj,mesh,l)
        %    gPar.type = obj.holeType;
        %    gPar.pnorm = obj.pnorm;   
        %    % a1 = obj.lattice.a1;
        %    % phi = atan2(a1(2), a1(1));
        %    % gPar.rotation = phi;
        %     coord = mesh.coord;
        %     xmin = min(coord(:,1));
        %     xmax = max(coord(:,1));
        %     ymin = min(coord(:,2));
        %     ymax = max(coord(:,2));          
        %     center_x = (xmin + xmax)/2;
        %     center_y = (ymin + ymax)/2;
        %     gPar.xCoorCenter = center_x;
        %     gPar.yCoorCenter = center_y;
        % 
        %     switch obj.meshType
        %         case 'Square'
        % 
        %              if size(coord,1) >= 4
        %                 v1 = coord(1,:);
        %                 v2 = coord(2,:);
        %                 v3 = coord(3,:);          
        %                 side1 = v2 - v1;
        %                 phi = atan2(side1(2), side1(1));
        %             else
        %                 phi = 0;
        %             end
        %         case 'Hexagon'
        % 
        %             if size(coord,1) >= 6
        %                 v1 = coord(1,:);
        %                 v2 = coord(2,:);
        %                 side1 = v2 - v1;
        %                 phi = atan2(side1(2), side1(1));
        %             else
        %                 phi = 0;
        %             end           
        %     end
        %     gPar.rotation = phi;
        % 
        %     switch obj.holeType
        %         case 'Circle'
        %             gPar.radius = l/2;
        %         case 'Square'
        %             gPar.length = l;
        % 
        %         case 'SmoothRectangle'
        % 
        %             sx = l;
        %             sy = l/2;            
        % 
        % 
        %             gPar.xSide = sx;
        %             gPar.ySide = sy;
        %             gPar.pnorm = 16;
        %         case 'Ellipse'
        %             gPar.type = "SmoothRectangle";
        %             gPar.xSide  = l(1);
        %             gPar.ySide  = l(2);
        %             gPar.pnorm  = 2;  
        %         case 'SmoothHexagon'
        %             gPar.radius = l;
        %             gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
        %         case 'ReinforcedHoneycomb'
        %             gPar.theta  = 1-l;                          
        %             gPar.eps    = 1;                        
        %             gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];            
        %             gPar.radius = l;    
        %             gPar.rotation = phi;
        %     end  
        %     g                  = GeometricalFunction(gPar);
        %     phiFun             = g.computeLevelSetFunction(mesh);
        %     lsCircle           = phiFun.fValues;
        %     % if l(1) <= 1e-9 && gPar.theta == 1
        %     %     ls = ones(size(lsCircle));
        %     % else
        %     ls = -lsCircle; 
        %     % end            
        % end

        function mat = createDensityMaterial(obj,lsf)
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1e-6*obj.E,obj.nu,obj.baseMesh.ndim);
            s.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1e-6*obj.E,obj.nu);
            s.matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(obj.E,obj.nu,obj.baseMesh.ndim);
            s.matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(obj.E,obj.nu);
            mI = MaterialInterpolator.create(s);

            x{1} = lsf;
            s.mesh                 = obj.baseMesh;
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = mI;
            s.dim                  = '2D';
            mat = Material.create(s);
        end

        function matHomog = solveElasticMicroProblem(obj,material,dens)           
            if obj.monitoring == true
                close all
                dens.plot
                shading interp
                colormap (flipud(pink))
                drawnow
            end

            s.mesh = obj.baseMesh;
            s.material = material;
            s.scale = 'MICRO';
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(obj.baseMesh);
            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            material.setDesignVariable({dens})
            fem.updateMaterial(material.obtainTensor())
            fem.solve();
            % fem.plotMicroFields(2,0);

            totVol = obj.baseMesh.computeVolume();
            matHomog = fem.Chomog/totVol;
        end

        function bc = createBoundaryConditions(obj,mesh)
            switch obj.meshType
                case 'Square'
                    
                    isCorner = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12) & ...
                                        (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
                    
                    sDir{1}.domain    = @(coor) isCorner(coor);
                    sDir{1}.direction = [1,2];
                    sDir{1}.value     = 0;
                    
                case 'Hexagon'
                   
                    isBottomVertex = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12) & ...
                                             (abs(coor(:,1) - min(coor(:,1))) < 1e-12);
                    
                    sDir{1}.domain    = @(coor) isBottomVertex(coor);
                    sDir{1}.direction = [1,2];
                    sDir{1}.value     = 0;
            end
        
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = 1;
            s.mesh = mesh;
            bc = BoundaryConditions(s);
            bc.updatePeriodicConditions(obj.masterSlave);
        end
        

       
        function fracVol = computeVolumeFraction(obj, b_val)
           
            rho = obj.createDensityLevelSet(b_val);    
            volDom = Integrator.compute(ConstantFunction.create(1, obj.baseMesh), obj.baseMesh, 2);
            fracVol = Integrator.compute(rho, rho.mesh, 2) / volDom;
        end
        
        function [mat] = assembleResults(obj,vec)
            sizeRes = size(vec);
            mat = zeros([sizeRes(1:end-1),obj.nSteps]);
            nStepsLastParam = obj.nSteps(end);
            nCombs = sizeRes(end);
            idxVec = repmat({':'}, 1, ndims(vec));
            idxMat = repmat({':'}, 1, ndims(mat));
            for i=1:nStepsLastParam
                idxVec{end} = i:nStepsLastParam:nCombs;
                idxMat{end} = i;
                mat(idxMat{:}) = vec(idxVec{:});
            end
        end
        
        function printTensorAtVolume(obj, targetVol)
        
            
            [~, idx] = min(abs(obj.volFrac - targetVol));
        
            fprintf('\nVolume solicitado: %.4f\n', targetVol);
            fprintf('Volume encontrado:  %.4f\n\n', obj.volFrac(idx));
        
            C = obj.Chomog(:,:,:,:,idx);
        
            
            Cvoigt = zeros(3,3);
        
            Cvoigt(1,1) = C(1,1,1,1);
            Cvoigt(1,2) = C(1,1,2,2);
            Cvoigt(1,3) = C(1,1,1,2);
        
            Cvoigt(2,1) = C(2,2,1,1);
            Cvoigt(2,2) = C(2,2,2,2);
            Cvoigt(2,3) = C(2,2,1,2);
        
            Cvoigt(3,1) = C(1,2,1,1);
            Cvoigt(3,2) = C(1,2,2,2);
            Cvoigt(3,3) = C(1,2,1,2);
        
            disp('Tensor homogenizado (Voigt):')
            disp(Cvoigt)
        
            
            K = (Cvoigt(1,1) + Cvoigt(2,2) + 2*Cvoigt(1,2))/4;
            mu = Cvoigt(3,3);
        
            fprintf('\nBulk modulus efetivo: %.6f\n', K);
            fprintf('Shear modulus efetivo: %.6f\n', mu);
        
           
            Eeff = mu*(3*K + mu)/(K + mu);
            nueff = (K - mu)/(K + mu);
        
            fprintf('Young efetivo: %.6f\n', Eeff);
            fprintf('Poisson efetivo: %.6f\n', nueff);
        
            
            figure;
            imagesc(Cvoigt)
            colorbar
            axis equal tight
            title(['Tensor homogenizado - vol = ', num2str(obj.volFrac(idx))])
        end



        function plot(obj)
            bMin = obj.paramHole(1);    
            bMax = obj.paramHole(end); 
            
            tiledlayout(1,3)
            
            nexttile
            hold on
            plot(obj.paramHole, squeeze(obj.Chomog(1,1,1,1,:)), 'LineStyle','none', 'Marker','o')
            fplot(obj.f{1,1,1,1}, [bMin, bMax])
            xlabel('b'); title('C_{1111}'); grid on
            
            nexttile
            hold on
            plot(obj.paramHole, squeeze(obj.Chomog(2,2,2,2,:)), 'LineStyle','none', 'Marker','o')
            fplot(obj.f{2,2,2,2}, [bMin, bMax])
            xlabel('b'); title('C_{2222}'); grid on
            
            nexttile
            hold on
            plot(obj.paramHole, squeeze(obj.Chomog(1,2,1,2,:)), 'LineStyle','none', 'Marker','o')
            fplot(obj.f{1,2,1,2}, [bMin, bMax])
            xlabel('b'); title('C_{1212}'); grid on
        end

        function fitting(obj)
            [obj.f,obj.df,obj.ddf] = DamageHomogenizationFitter.computePolynomial(12,obj.paramHole,obj.Chomog);
        end
        
    end

   
    
end

