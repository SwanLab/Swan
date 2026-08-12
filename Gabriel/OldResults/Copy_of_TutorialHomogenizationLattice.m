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

        %% ========== MÉTODOS DE VALIDAÇÃO ==========

        function validateHomogenization(obj, b_values)
            % Valida os resultados da homogeneização usando o fitting existente
            % b_values: vetor de valores de b para testar (ex: [0, 0.2, 0.4, 0.6, 0.8, 1.0])
            
            if nargin < 2
                b_values = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
            end
            
            nTests = length(b_values);
            
            fprintf('\n========== VALIDAÇÃO DA HOMOGENEIZAÇÃO ==========\n\n');
            
            % 1. TABELA DE RESULTADOS
            fprintf('Tabela 1: Componentes do tensor homogeneizado\n');
            fprintf('%-10s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
                    'b', 'C11', 'C22', 'C33', 'C12', 'C11/C22', 'Zener');
            fprintf('%s\n', repmat('-', 1, 85));
            
            for i = 1:nTests
                b_val = b_values(i);
                % Encontrar índice mais próximo
                [~, idx] = min(abs(obj.paramHole - b_val));
                
                C11 = obj.Chomog(1,1,1,1,idx);
                C22 = obj.Chomog(2,2,2,2,idx);
                C33 = obj.Chomog(1,2,1,2,idx);
                C12 = obj.Chomog(1,1,2,2,idx);
                Zener = 2*C33/(C11 - C12);
                
                fprintf('%-10.4f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', ...
                        b_val, C11, C22, C33, C12, C11/C22, Zener);
            end
            
            % 2. DERIVADAS USANDO O FITTING
            if ~isempty(obj.df)
                fprintf('\nTabela 2: Derivadas dos componentes (dC/db) via fitting\n');
                fprintf('%-10s %-12s %-12s %-12s\n', 'b', 'dC11/db', 'dC22/db', 'dC33/db');
                fprintf('%s\n', repmat('-', 1, 50));
                
                for i = 1:nTests
                    b_val = b_values(i);
                    
                    % Usar as funções derivadas do fitting
                    dC11 = obj.df{1,1,1,1}(b_val);
                    dC22 = obj.df{2,2,2,2}(b_val);
                    dC33 = obj.df{1,2,1,2}(b_val);
                    
                    fprintf('%-10.4f %-12.6f %-12.6f %-12.6f\n', b_val, dC11, dC22, dC33);
                end
            end
            
            fprintf('\n================================================\n');
        end

        function plotDerivativesComparison(obj)
            % Compara derivadas via diferenças finitas vs derivadas do fitting
            
            b_vals = obj.paramHole;
            n = length(b_vals);
            
            % Derivadas por diferenças finitas
            dC11_fd = zeros(n,1);
            dC22_fd = zeros(n,1);
            dC33_fd = zeros(n,1);
            
            for i = 2:n-1
                db = b_vals(i+1) - b_vals(i-1);
                dC11_fd(i) = (obj.Chomog(1,1,1,1,i+1) - obj.Chomog(1,1,1,1,i-1)) / db;
                dC22_fd(i) = (obj.Chomog(2,2,2,2,i+1) - obj.Chomog(2,2,2,2,i-1)) / db;
                dC33_fd(i) = (obj.Chomog(1,2,1,2,i+1) - obj.Chomog(1,2,1,2,i-1)) / db;
            end
            
            figure('Name', 'Comparação de Derivadas', 'Position', [100, 100, 1200, 400]);
            
            subplot(1,3,1);
            plot(b_vals, dC11_fd, 'bo-', 'LineWidth', 1.5); hold on;
            fplot(@(x) obj.df{1,1,1,1}(x), [b_vals(1), b_vals(end)], 'r-', 'LineWidth', 1.5);
            xlabel('b'); ylabel('dC_{11}/db'); title('Derivada C11');
            legend('Diferenças finitas', 'Fitting', 'Location', 'best');
            grid on;
            
            subplot(1,3,2);
            plot(b_vals, dC22_fd, 'bo-', 'LineWidth', 1.5); hold on;
            fplot(@(x) obj.df{2,2,2,2}(x), [b_vals(1), b_vals(end)], 'r-', 'LineWidth', 1.5);
            xlabel('b'); ylabel('dC_{22}/db'); title('Derivada C22');
            legend('Diferenças finitas', 'Fitting', 'Location', 'best');
            grid on;
            
            subplot(1,3,3);
            plot(b_vals, dC33_fd, 'bo-', 'LineWidth', 1.5); hold on;
            fplot(@(x) obj.df{1,2,1,2}(x), [b_vals(1), b_vals(end)], 'r-', 'LineWidth', 1.5);
            xlabel('b'); ylabel('dC_{33}/db'); title('Derivada C33');
            legend('Diferenças finitas', 'Fitting', 'Location', 'best');
            grid on;
            
            sgtitle('Comparação de Derivadas');
        end

        function createValidationTable(obj)
            % Cria uma tabela completa de validação com todos os b
            
            b_vals = obj.paramHole;
            n = length(b_vals);
            
            results = zeros(n, 7);  % b, C11, C22, C33, C12, dC11, dC22
            
            for i = 1:n
                results(i,1) = b_vals(i);
                results(i,2) = obj.Chomog(1,1,1,1,i);
                results(i,3) = obj.Chomog(2,2,2,2,i);
                results(i,4) = obj.Chomog(1,2,1,2,i);
                results(i,5) = obj.Chomog(1,1,2,2,i);
                
                if ~isempty(obj.df)
                    results(i,6) = obj.df{1,1,1,1}(b_vals(i));
                    results(i,7) = obj.df{2,2,2,2}(b_vals(i));
                end
            end
            
            % Salvar tabela
            T = array2table(results, 'VariableNames', {'b', 'C11', 'C22', 'C33', 'C12', 'dC11', 'dC22'});
            writetable(T, 'validation_results.csv');
            disp(T);
            
            % Plotar
            figure('Name', 'Validação Completa', 'Position', [50, 50, 1400, 800]);
            
            subplot(2,2,1);
            plot(results(:,1), results(:,2), 'b-o', results(:,1), results(:,3), 'r-o', results(:,1), results(:,4), 'g-o');
            xlabel('b'); ylabel('C'); legend('C11', 'C22', 'C33', 'Location', 'best'); title('Componentes do tensor'); grid on;
            
            if ~isempty(obj.df)
                subplot(2,2,2);
                plot(results(:,1), results(:,6), 'b-o', results(:,1), results(:,7), 'r-o');
                xlabel('b'); ylabel('dC/db'); legend('dC11', 'dC22', 'Location', 'best'); title('Derivadas (fitting)'); grid on;
            end
            
            subplot(2,2,3);
            plot(results(:,1), results(:,2)./results(:,3), 'k-o');
            xlabel('b'); ylabel('C11/C22'); title('Anisotropia'); grid on; yline(1, 'r--', 'Isotropia');
            
            subplot(2,2,4);
            Zener = 2*results(:,4) ./ (results(:,2) - results(:,5));
            plot(results(:,1), Zener, 'm-o');
            xlabel('b'); ylabel('Zener Ratio'); title('Anisotropia de cisalhamento'); grid on; yline(1, 'r--', 'Isotropia');
            
            sgtitle('Validação dos Resultados de Homogeneização');
        end

        function plotFieldsForB(obj, b_val, scale)
            % Plota campos de deslocamento e deformação para um b específico
            % b_val: valor de b
            % scale: fator de escala para deformação (padrão=50)
            
            if nargin < 3
                scale = 50;
            end
            
            % Encontrar índice
            [~, idx] = min(abs(obj.paramHole - b_val));
            b_actual = obj.paramHole(idx);
            
            fprintf('\nPlotando campos para b = %.4f\n', b_actual);
            
            % Recalcular o problema para este b
            a_val = exp(b_actual^2) - 0.5;
            d_val = (1 + b_actual^2) / a_val;
            
            v1 = [a_val, b_actual];
            v2 = [b_actual, d_val];
            obj.latticeVectors = [v1; v2];
            obj.defineMesh();
            
            dens = obj.createDensityLevelSet(obj.fixedHoleSize);
            mat = obj.createDensityMaterial(dens);
            
            s.mesh = obj.baseMesh;
            s.material = mat;
            s.scale = 'MICRO';
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(obj.baseMesh);
            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            mat.setDesignVariable({dens})
            fem.updateMaterial(mat.obtainTensor())
            fem.solve();
            
            % Plotar para cada base
            figure('Name', sprintf('Campos para b = %.4f', b_actual), 'Position', [100, 100, 1400, 1000]);
            
            for iBasis = 1:3
                % Deslocamentos flutuantes
                uFlucFun = fem.uFluc{iBasis};
                uFluc = uFlucFun.fValues;
                
                switch iBasis
                    case 1
                        E_macro = [1, 0; 0, 0];
                        title_str = 'Base 1: Tração Horizontal (ε_{xx}=1)';
                    case 2
                        E_macro = [0, 0; 0, 1];
                        title_str = 'Base 2: Tração Vertical (ε_{yy}=1)';
                    case 3
                        E_macro = [0, 0.5; 0.5, 0];
                        title_str = 'Base 3: Cisalhamento (ε_{xy}=1)';
                end
                
                coord = obj.baseMesh.coord;
                u_macro = (E_macro * coord')';
                u_total = u_macro + uFluc;
                coord_def = coord + scale * u_total;
                
                % Subplot 1: Deslocamento X
                subplot(3,3, (iBasis-1)*3 + 1);
                uPlot = LagrangianFunction.create(obj.baseMesh, 1, 'P1');
                uPlot.setFValues(uFluc(:,1));
                uPlot.plot();
                title(sprintf('%s - u_x flutuante', title_str));
                colorbar;
                axis equal;
                
                % Subplot 2: Deslocamento Y
                subplot(3,3, (iBasis-1)*3 + 2);
                uPlot = LagrangianFunction.create(obj.baseMesh, 1, 'P1');
                uPlot.setFValues(uFluc(:,2));
                uPlot.plot();
                title(sprintf('%s - u_y flutuante', title_str));
                colorbar;
                axis equal;
                
                % Subplot 3: Malha deformada
                subplot(3,3, (iBasis-1)*3 + 3);
                trisurf(obj.baseMesh.connec, coord_def(:,1), coord_def(:,2), zeros(size(coord_def,1),1), ...
                        'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5);
                view(2);
                axis equal;
                title(sprintf('%s - Malha deformada (escala=%d)', title_str, scale));
                xlabel('x'); ylabel('y');
            end
            
            sgtitle(sprintf('Campos da célula periódica - b = %.4f, a = %.4f, d = %.4f', b_actual, a_val, d_val));
        end

        function checkPeriodicCondition(obj, b_val)
            % Verifica se as condições periódicas estão sendo satisfeitas
            % Compara deslocamentos flutuantes em pares master-slave
            
            [~, idx] = min(abs(obj.paramHole - b_val));
            b_actual = obj.paramHole(idx);
            
            % Recalcular o problema
            a_val = exp(b_actual^2) - 0.5;
            d_val = (1 + b_actual^2) / a_val;
            
            v1 = [a_val, b_actual];
            v2 = [b_actual, d_val];
            obj.latticeVectors = [v1; v2];
            obj.defineMesh();
            
            dens = obj.createDensityLevelSet(obj.fixedHoleSize);
            mat = obj.createDensityMaterial(dens);
            
            s.mesh = obj.baseMesh;
            s.material = mat;
            s.scale = 'MICRO';
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(obj.baseMesh);
            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            mat.setDesignVariable({dens})
            fem.updateMaterial(mat.obtainTensor())
            fem.solve();
            
            % Verificar periodicidade
            masterNodes = obj.masterSlave(:,1);
            slaveNodes = obj.masterSlave(:,2);
            
            figure('Name', sprintf('Verificação de Periodicidade - b = %.4f', b_actual), 'Position', [100, 100, 1400, 500]);
            
            for iBasis = 1:3
                uFlucFun = fem.uFluc{iBasis};
                uFluc = uFlucFun.fValues;
                
                u_master = uFluc(masterNodes, :);
                u_slave = uFluc(slaveNodes, :);
                
                diff_u = sqrt((u_slave(:,1) - u_master(:,1)).^2 + (u_slave(:,2) - u_master(:,2)).^2);
                
                subplot(1,3,iBasis);
                semilogy(diff_u, 'o-', 'MarkerSize', 4);
                xlabel('Par master-slave');
                ylabel('||u_{slave} - u_{master}||');
                title(sprintf('Base %d - Erro periódico médio: %.2e', iBasis, mean(diff_u)));
                grid on;
            end
            
            sgtitle(sprintf('Verificação de Periodicidade - b = %.4f (u^{fluc}_{slave} - u^{fluc}_{master} = 0)', b_actual));
        end

        function runValidationSuite(obj)
            % Executa toda a suíte de validação
            
            fprintf('\n');
            fprintf('╔════════════════════════════════════════════════════════════════╗\n');
            fprintf('║                    SUÍTE DE VALIDAÇÃO                           ║\n');
            fprintf('╚════════════════════════════════════════════════════════════════╝\n');
            
            % 1. Tabela de resultados
            obj.validateHomogenization([0, 0.2, 0.4, 0.6, 0.8, 1.0]);
            
            % 2. Comparação de derivadas
            if ~isempty(obj.df)
                obj.plotDerivativesComparison();
            end
            
            % 3. Tabela completa
            obj.createValidationTable();
            
            % 4. Verificação de periodicidade para alguns valores
            fprintf('\nVerificando periodicidade para b = 0, 0.5, 1.0...\n');
            obj.checkPeriodicCondition(0);
            obj.checkPeriodicCondition(0.5);
            obj.checkPeriodicCondition(1.0);
            
            fprintf('\n✅ Validação concluída!\n');
        end

    end
    
    methods (Access = private)
        
        function init(obj)
            obj.E          = 1;
            obj.nu         = 0.3;
            obj.meshType   = 'Square';
            obj.meshN      = 100;

            obj.holeType   = 'Square';
            obj.pnorm      = 'Inf';
            % obj.damageType = 'Area';
            obj.nSteps     = 60;
            obj.monitoring = true;
            obj.fixedHoleSize = 0.7;
        end

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

        function computeHoleParams(obj)
            obj.maxParam = 1;
            obj.paramHole = linspace(-0.5,obj.maxParam, obj.nSteps);
        end

        function compute(obj)
            nComb = length(obj.paramHole);
            mat = zeros(2,2,2,2,nComb);
            volF = zeros(1,nComb);
            
            for i = 1:nComb
                b_val = obj.paramHole(i);
                obj.currentB = b_val;             
                
                a = 1 + 0.45*b_val + 0.20*b_val^2;
                d = (1 + b_val^2) / a;
                
                v1 = [a,    b_val];
                v2 = [b_val, d   ];
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
            [obj.f,obj.df,obj.ddf] = DamageHomogenizationFitter.computePolynomial(9,obj.paramHole,obj.Chomog);
        end
        
    end

end