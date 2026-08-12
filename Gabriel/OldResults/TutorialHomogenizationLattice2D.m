classdef TutorialHomogenizationLattice2D < handle
    % Tutorial para explorar a matriz simétrica T = [a, b; b, (1+b²)/a]
    
    properties (Access = public)
        Chomog
        volFrac
    end

    properties (Access = private)
        E
        nu
        meshType
        meshN
        holeType
        pnorm
        monitoring
        Mmass
        fixedHoleSize
        latticeVectors
        baseMesh
        masterSlave
        test
    end

    methods (Access = public)
        
        function obj = TutorialHomogenizationLattice2D()
            obj.init();
        end
        function exportSurfaceData(obj, all_a, all_b, all_C11, all_C22, all_C33, all_C12, vf_target, filename)
            if nargin < 9
                filename = sprintf('Chomog_Vf%.2f.csv', vf_target);
            end
            
            [nA, nB] = size(all_C11);
            nCases   = nA * nB;
            
            % Criar matriz de dados
            data = zeros(nCases, 8);
            idx  = 0;
            
            for iA = 1:nA
                for iB = 1:nB
                    idx = idx + 1;
                    a   = all_a(iA, iB);
                    b   = all_b(iA, iB);
                    d   = (1 + b^2) / a;
                    
                    data(idx, 1) = a;
                    data(idx, 2) = b;
                    data(idx, 3) = d;
                    data(idx, 4) = all_C11(iA, iB);
                    data(idx, 5) = all_C22(iA, iB);
                    data(idx, 6) = all_C33(iA, iB);
                    data(idx, 7) = all_C12(iA, iB);
                    data(idx, 8) = all_C11(iA,iB) / all_C22(iA,iB);  % anisotropia
                end
            end
            
            % Escrever CSV com header
            fid = fopen(filename, 'w');
            fprintf(fid, 'a,b,d,C11,C22,C33,C12,anisotropy\n');
            for i = 1:nCases
                fprintf(fid, '%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', data(i,:));
            end
            fclose(fid);
            
            fprintf('Exportado: %s (%d casos)\n', filename, nCases);
        end
        
        function exploreSymmetricCell(obj, a_vals, b_vals)
            % Explora células com T = [a, b; b, (1+b²)/a]
            
            if nargin < 2
                a_vals = [0.6, 0.8, 1.0, 1.2, 1.4];
                b_vals = [-1.0, -0.5, 0, 0.5, 1.0];
            end
            
            nA = length(a_vals);
            nB = length(b_vals);
            
            fprintf('\n========================================\n');
            fprintf('Explorando matriz simétrica\n');
            fprintf('T = [a, b; b, (1+b²)/a]\n');
            fprintf('det(T) = 1 (área constante)\n');
            fprintf('========================================\n\n');
            
            % Armazenar resultados
            all_a = zeros(nA, nB);
            all_b = zeros(nA, nB);
            all_C11 = zeros(nA, nB);
            all_C22 = zeros(nA, nB);
            all_C33 = zeros(nA, nB);
            all_C12 = zeros(nA, nB);
            all_meshes = cell(nA, nB);
            
            for iA = 1:nA
                for iB = 1:nB
                    a = a_vals(iA);
                    b = b_vals(iB);
                    d = (1 + b^2) / a;
                    
                    fprintf('a=%.2f, b=%.2f, d=%.4f\n', a, b, d);
                    
                    % Construir latticeVectors
                    v1 = [a, b];
                    v2 = [b, d];
                    obj.latticeVectors = [v1; v2];
                    obj.defineMesh();
                    
                    % Resolver homogeneização (já normalizado)
                    C = obj.computeHomogenization();
                    
                    % Guardar
                    all_a(iA, iB) = a;
                    all_b(iA, iB) = b;
                    all_C11(iA, iB) = C(1,1,1,1);
                    all_C22(iA, iB) = C(2,2,2,2);
                    all_C33(iA, iB) = C(1,2,1,2);
                    all_C12(iA, iB) = C(1,1,2,2);
                    all_meshes{iA, iB} = obj.baseMesh;
                end
            end
            
            % Plotar células
            obj.plotCells(all_meshes, all_a, all_b);
            
            % Plotar resultados
            obj.plotResults(all_a, all_b, all_C11, all_C22, all_C33, all_C12);
        end
        function plotSurfacesForVf(obj, vf_target)
            if nargin < 2
                vf_target = 0.3;
            end
            
            % Calcular holeSize para o Vf desejado
            obj.fixedHoleSize = sqrt(1 - vf_target);
            fprintf('Vf = %.2f  →  holeSize = %.4f\n', vf_target, obj.fixedHoleSize);
            
            % Grelha fina de (a,b)
            a_vals = linspace(0.7, 2, 30);
            b_vals = linspace(-1, 1, 30);
            
            nA = length(a_vals);
            nB = length(b_vals);
            
            all_C11 = zeros(nA, nB);
            all_C22 = zeros(nA, nB);
            all_C33 = zeros(nA, nB);
            all_C12 = zeros(nA, nB);
            all_meshes = cell(nA, nB);
            all_a = zeros(nA, nB);
            all_b = zeros(nA, nB);
            
            nTotal = nA * nB;
            idx = 0;
            
            for iA = 1:nA
                for iB = 1:nB
                    idx = idx + 1;
                    a = a_vals(iA);
                    b = b_vals(iB);
                    d = (1 + b^2) / a;
                    
                    fprintf('  %d/%d — a=%.2f, b=%.2f\n', idx, nTotal, a, b);
                    
                    v1 = [a, b];
                    v2 = [b, d];
                    obj.latticeVectors = [v1; v2];
                    obj.defineMesh();
                    
                    C = obj.computeHomogenization();
                    
                    all_a(iA, iB)   = a;
                    all_b(iA, iB)   = b;
                    all_C11(iA, iB) = C(1,1,1,1);
                    all_C22(iA, iB) = C(2,2,2,2);
                    all_C33(iA, iB) = C(1,2,1,2);
                    all_C12(iA, iB) = C(1,1,2,2);
                    all_meshes{iA, iB} = obj.baseMesh;
                end
            end
            obj.exportSurfaceData(all_a, all_b, all_C11, all_C22, all_C33, all_C12, vf_target);
            % Plot das células com furos
            obj.plotCellsWithHoles(all_meshes, all_a, all_b, vf_target);
            
            % Plot das superfícies
            obj.plotSurfaces(all_a, all_b, all_C11, all_C22, all_C33, all_C12, vf_target);
        end
        function plotTensorVsB(obj, vf_target, a_vals, b_vals)
            if nargin < 2; vf_target = 0.9; end
            if nargin < 3; a_vals = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8]; end
            if nargin < 4; b_vals = linspace(-1, 1, 30); end
            
            obj.fixedHoleSize = sqrt(1 - vf_target);
            fprintf('Vf = %.2f  →  holeSize = %.4f\n', vf_target, obj.fixedHoleSize);
            
            nA   = length(a_vals);
            nB   = length(b_vals);
            colors = lines(nA);
            
            % Armazenar resultados
            C11 = zeros(nA, nB);
            C22 = zeros(nA, nB);
            C33 = zeros(nA, nB);
            C12 = zeros(nA, nB);
            
            nTotal = nA * nB;
            idx = 0;
            
            for iA = 1:nA
                for iB = 1:nB
                    idx = idx + 1;
                    a = a_vals(iA);
                    b = b_vals(iB);
                    d = (1 + b^2) / a;
                    
                    fprintf('  %d/%d — a=%.2f, b=%.2f\n', idx, nTotal, a, b);
                    
                    v1 = [a, b];
                    v2 = [b, d];
                    obj.latticeVectors = [v1; v2];
                    obj.defineMesh();
                    
                    C = obj.computeHomogenization();
                    
                    C11(iA, iB) = C(1,1,1,1);
                    C22(iA, iB) = C(2,2,2,2);
                    C33(iA, iB) = C(1,2,1,2);
                    C12(iA, iB) = C(1,1,2,2);
                end
                
            end
            
            % Mascarar zeros (casos inválidos)
            C11(C11 < 1e-3) = NaN;
            C22(C22 < 1e-3) = NaN;
            C33(C33 < 1e-3) = NaN;
            C12(C12 < 1e-3) = NaN;
            obj.plotMicrosForCurves(a_vals, b_vals, vf_target);
            % Plot
            figure('Name', sprintf('C^h vs b — Vf=%.2f', vf_target), ...
                   'Position', [50, 50, 1400, 900]);
            
            labels = arrayfun(@(a) sprintf('a=%.1f', a), a_vals, 'UniformOutput', false);
            
            subplot(2,3,1); hold on;
            for iA = 1:nA
                plot(b_vals, C11(iA,:), '-', 'Color', colors(iA,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('C_{11}'); title('C_{11} vs b');
            legend(labels, 'Location', 'best'); grid on;
            
            subplot(2,3,2); hold on;
            for iA = 1:nA
                plot(b_vals, C22(iA,:), '-', 'Color', colors(iA,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('C_{22}'); title('C_{22} vs b');
            legend(labels, 'Location', 'best'); grid on;
            
            subplot(2,3,3); hold on;
            for iA = 1:nA
                plot(b_vals, C33(iA,:), '-', 'Color', colors(iA,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('C_{33}'); title('C_{33} vs b');
            legend(labels, 'Location', 'best'); grid on;
            
            subplot(2,3,4); hold on;
            for iA = 1:nA
                plot(b_vals, C12(iA,:), '-', 'Color', colors(iA,:), 'LineWidth', 1.5);
            end
            xlabel('b'); ylabel('C_{12}'); title('C_{12} vs b');
            legend(labels, 'Location', 'best'); grid on;
            
            subplot(2,3,5); hold on;
            for iA = 1:nA
                plot(b_vals, C11(iA,:)./C22(iA,:), '-', 'Color', colors(iA,:), 'LineWidth', 1.5);
            end
            yline(1, 'k--', 'Isotropia');
            xlabel('b'); ylabel('C_{11}/C_{22}'); title('Anisotropia vs b');
            legend(labels, 'Location', 'best'); grid on;
            
            subplot(2,3,6); hold on;
            for iA = 1:nA
                denom = C11(iA,:) - C12(iA,:);
                denom(abs(denom) < 1e-10) = NaN;
                plot(b_vals, 2*C33(iA,:)./denom, '-', 'Color', colors(iA,:), 'LineWidth', 1.5);
            end
            yline(1, 'k--', 'Isotropia');
            xlabel('b'); ylabel('A'); title('Zener Ratio vs b');
            legend(labels, 'Location', 'best'); grid on;
            
            sgtitle(sprintf('C^h(a,b) — Vf=%.2f, curvas para diferentes a', vf_target));
        end
        function exploreCurveAofB(obj, vf_target, b_vals)
            if nargin < 2; vf_target = 0.6;               end
            if nargin < 3; b_vals = linspace(-1, 1, 30);  end
            
            obj.fixedHoleSize = sqrt(1 - vf_target);
            fprintf('Vf = %.2f  →  holeSize = %.4f\n', vf_target, obj.fixedHoleSize);
            
            nB   = length(b_vals);
            C11  = zeros(1, nB);
            C22  = zeros(1, nB);
            C33  = zeros(1, nB);
            C12  = zeros(1, nB);
            a_computed = zeros(1, nB);
            
            for iB = 1:nB
                b = b_vals(iB);
                a = (exp(b^2) - 1) + 1 ;      % <-- relação a(b)
                d = (1 + b^2) / a;      % = sqrt(1+b²) = a  →  T diagonal igual!
                
                a_computed(iB) = a;
                fprintf('  b=%.2f → a=%.4f, d=%.4f\n', b, a, d);
                
                obj.latticeVectors = [a, b; b, d];
                obj.defineMesh();
                C = obj.computeHomogenization();
                
                C11(iB) = C(1,1,1,1);
                C22(iB) = C(2,2,2,2);
                C33(iB) = C(1,2,1,2);
                C12(iB) = C(1,1,2,2);
            end
            
            % Mascarar zeros
            C11(C11 < 1e-3) = NaN;
            C22(C22 < 1e-3) = NaN;
            C33(C33 < 1e-3) = NaN;
            C12(C12 < 1e-3) = NaN;
            
            % Plot das micros primeiro
            obj.plotMicrosOneCurve(b_vals, a_computed, vf_target);
            
            % Plot das curvas
            figure('Name', sprintf('a(b)=sqrt(1+b²) — Vf=%.2f', vf_target), ...
                   'Position', [50, 50, 1400, 800]);
            
            subplot(2,3,1);
            plot(b_vals, a_computed, 'k-', 'LineWidth', 2);
            xlabel('b'); ylabel('a(b)'); title('Relação a(b) = \sqrt{1+b^2}'); grid on;
            
            subplot(2,3,2); hold on;
            plot(b_vals, C11, 'b-', 'LineWidth', 2, 'DisplayName', 'C_{11}');
            plot(b_vals, C22, 'r-', 'LineWidth', 2, 'DisplayName', 'C_{22}');
            plot(b_vals, C33, 'g-', 'LineWidth', 2, 'DisplayName', 'C_{33}');
            legend; xlabel('b'); title('Componentes vs b'); grid on;
            
            subplot(2,3,3);
            plot(b_vals, C12, 'k-', 'LineWidth', 2);
            xlabel('b'); ylabel('C_{12}'); title('C_{12} vs b'); grid on;
            
            subplot(2,3,4);
            plot(b_vals, C11./C22, 'k-', 'LineWidth', 2);
            yline(1, 'r--', 'Isotropia');
            xlabel('b'); ylabel('C_{11}/C_{22}'); title('Anisotropia vs b'); grid on;
            
            subplot(2,3,5);
            denom = C11 - C12;
            denom(abs(denom) < 1e-10) = NaN;
            plot(b_vals, 2*C33./denom, 'k-', 'LineWidth', 2);
            yline(1, 'r--', 'Isotropia');
            xlabel('b'); ylabel('A'); title('Zener Ratio vs b'); grid on;
            
            subplot(2,3,6);
            axis off;
            text(0.1, 0.9, sprintf('Vf = %.2f', vf_target), 'FontSize', 12);
            text(0.1, 0.8, 'a(b) = \sqrt{1+b^2}', 'FontSize', 11);
            text(0.1, 0.7, 'T = [a,b; b,a]  (diagonal igual!)', 'FontSize', 11);
            text(0.1, 0.6, 'det(T) = a^2 - b^2 = 1', 'FontSize', 11);
            
            sgtitle(sprintf('C^h com a(b) = sqrt(1+b²) — Vf=%.2f', vf_target));
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.E          = 1;
            obj.nu         = 0.3;
            obj.meshType   = 'Square';
            obj.meshN      = 70;
            obj.holeType   = 'Square';
            obj.pnorm      = 'Inf';
            obj.monitoring = false;
            obj.fixedHoleSize = 0.7;
        end
        
        function defineMesh(obj)
            s.latticeVectors = obj.latticeVectors;
            s.divUnit = obj.meshN;
            s.filename = '';
            MC = MeshCreator(s);
            MC.computeMeshNodes();
            
            s.coord  = MC.coord;
            s.connec = MC.connec;
            obj.baseMesh = Mesh.create(s);
            obj.masterSlave = MC.masterSlaveIndex;
            obj.test = LagrangianFunction.create(obj.baseMesh,1,'P1');
            obj.Mmass = IntegrateLHS(@(u,v) DP(v,u), obj.test, obj.test, obj.baseMesh, 'Domain');
        end
        
        function matHomog = computeHomogenization(obj)
            dens = obj.createDensityLevelSet();
            mat  = obj.createDensityMaterial(dens);
            matHomog = obj.solveElasticMicroProblem(mat, dens);
        end
        
        function lsf = createDensityLevelSet(obj)
            ls = obj.computeLevelSet(obj.baseMesh, obj.fixedHoleSize);
            sUm.backgroundMesh = obj.baseMesh;
            sUm.boundaryMesh   = obj.baseMesh.createBoundaryMesh;
            uMesh = UnfittedMesh(sUm);
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
                gPar.a1 = [1, 0];
                gPar.a2 = [0, 1];
                phi = 0;
            end
            gPar.rotation = phi;
            
            gPar.length = l_fixo;
            
            g = GeometricalFunction(gPar);
            phiFun = g.computeLevelSetFunction(mesh);
            lsCircle = phiFun.fValues;
            ls = -lsCircle;
        end
        
        function mat = createDensityMaterial(obj, lsf)
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1e-6*obj.E, obj.nu, obj.baseMesh.ndim);
            s.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1e-6*obj.E, obj.nu);
            s.matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(obj.E, obj.nu, obj.baseMesh.ndim);
            s.matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(obj.E, obj.nu);
            mI = MaterialInterpolator.create(s);
            
            x{1} = lsf;
            s.mesh                 = obj.baseMesh;
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = mI;
            s.dim                  = '2D';
            mat = Material.create(s);
        end
        
        function matHomog = solveElasticMicroProblem(obj, material, dens)
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
            
            totVol = obj.baseMesh.computeVolume();
            matHomog = fem.Chomog / totVol;
        end
        
        function bc = createBoundaryConditions(obj, mesh)
            isCorner = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12) & ...
                                (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
            
            sDir{1}.domain    = @(coor) isCorner(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
            
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
        
        function plotCells(obj, all_meshes, all_a, all_b)
            [nA, nB] = size(all_meshes);
            
            % Selecionar índices para plotar (9 casos representativos)
            idxA = unique([1, round(nA/2), nA]);
            idxB = unique([1, round(nB/2), nB]);
            
            nPlotA = length(idxA);
            nPlotB = length(idxB);
            
            figure('Name', 'Celdas unitárias', 'Position', [100, 100, 300*nPlotB, 300*nPlotA]);
            
            plotIdx = 0;
            for iA = idxA
                for iB = idxB
                    plotIdx = plotIdx + 1;
                    
                    mesh = all_meshes{iA, iB};
                    a = all_a(iA, iB);
                    b = all_b(iA, iB);
                    
                    subplot(nPlotA, nPlotB, plotIdx);
                    patch('Faces', mesh.connec, 'Vertices', mesh.coord, ...
                          'EdgeColor', [0.2 0.4 0.8], 'FaceColor', 'none', 'LineWidth', 0.3);
                    axis equal; axis off;
                    title(sprintf('a=%.2f, b=%.2f', a, b), 'FontSize', 9);
                    
                    % Desenhar vetores
                    hold on;
                    v1 = [a, b];
                    v2 = [b, (1+b^2)/a];
                    origin = [min(mesh.coord(:,1)), min(mesh.coord(:,2))];
                    scale = 0.25;
                    quiver(origin(1), origin(2), v1(1)*scale, v1(2)*scale, 0, ...
                           'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
                    quiver(origin(1), origin(2), v2(1)*scale, v2(2)*scale, 0, ...
                           'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
                    hold off;
                end
            end
            
            sgtitle('Celdas: T = [a, b; b, (1+b²)/a] — rojo: v1, verde: v2');
        end
        
        function plotResults(obj, all_a, all_b, C11, C22, C33, C12)
            
            % Obter os valores únicos de a e b (primeira linha/coluna)
            a_unique = all_a(:,1);  % vector de a
            b_unique = all_b(1,:);  % vector de b
            
            figure('Name', 'Tensor homogeneizado', 'Position', [50, 50, 1400, 800]);
            
            % C11
            subplot(2,3,1);
            surf(b_unique, a_unique, C11);
            xlabel('b'); ylabel('a'); zlabel('C_{11}');
            title('C_{11}'); colorbar; shading interp;
            
            % C22
            subplot(2,3,2);
            surf(b_unique, a_unique, C22);
            xlabel('b'); ylabel('a'); zlabel('C_{22}');
            title('C_{22}'); colorbar; shading interp;
            
            % C33
            subplot(2,3,3);
            surf(b_unique, a_unique, C33);
            xlabel('b'); ylabel('a'); zlabel('C_{33}');
            title('C_{33} (shear)'); colorbar; shading interp;
            
            % Anisotropia
            subplot(2,3,4);
            surf(b_unique, a_unique, C11./C22);
            xlabel('b'); ylabel('a'); zlabel('C_{11}/C_{22}');
            title('Anisotropia'); colorbar; shading interp;
            hold on;
            
            % Zener
            subplot(2,3,5);
            denom = C11 - C12;
            denom(abs(denom) < 1e-10) = NaN;
            surf(b_unique, a_unique, 2*C33./denom);
            xlabel('b'); ylabel('a'); zlabel('A');
            title('Zener Ratio'); colorbar; shading interp;
            
            % Info
            subplot(2,3,6);
            axis off;
            text(0.1, 0.9, sprintf('Vf = %.3f', 1-obj.fixedHoleSize^2), 'FontSize', 12);
            text(0.1, 0.8, 'T = [a,b; b,(1+b²)/a]', 'FontSize', 11);
            text(0.1, 0.7, 'det(T) = 1', 'FontSize', 11);
            
            sgtitle('C^h(a,b) — superfície de resposta');
        end
        function plotCellsWithHoles(obj, all_meshes, all_a, all_b, vf_target)
            [nA, nB] = size(all_meshes);
            
            % 9 casos representativos
            idxA = unique([1, round(nA/2), nA]);
            idxB = unique([1, round(nB/2), nB]);
            nPlotA = length(idxA);
            nPlotB = length(idxB);
            
            figure('Name', sprintf('Células com furos — Vf=%.2f', vf_target), ...
                   'Position', [100, 100, 350*nPlotB, 350*nPlotA]);
            
            plotIdx = 0;
            for iA = idxA
                for iB = idxB
                    plotIdx = plotIdx + 1;
                    
                    % Usar a malha já guardada
                    mesh = all_meshes{iA, iB};
                    a    = all_a(iA, iB);
                    b    = all_b(iA, iB);
                    
                    % Recalcular só a densidade — sem redefinir a malha global
                    obj.baseMesh      = mesh;
                    obj.latticeVectors = [a, b; b, (1+b^2)/a];
                    obj.test          = LagrangianFunction.create(mesh, 1, 'P1');
                    dens              = obj.createDensityLevelSet();
                    
                    subplot(nPlotA, nPlotB, plotIdx);
                    
                    % Plotar manualmente sem abrir nova figura
                    rho    = dens.project('P0');
                    rhoVal = squeeze(rho.fValues);
                    patch('Faces', mesh.connec, 'Vertices', mesh.coord, ...
                          'FaceVertexCData', rhoVal, 'FaceColor', 'flat', ...
                          'EdgeColor', 'none');
                    colormap(flipud(gray));
                    caxis([0 1]);
                    axis equal; axis off;
                    title(sprintf('a=%.2f, b=%.2f', a, b), 'FontSize', 9);
                end
            end
            
            sgtitle(sprintf('Células com furos — Vf=%.2f', vf_target));
        end

        function plotSurfaces(obj, all_a, all_b, C11, C22, C33, C12, vf_target)
            a_unique = all_a(:,1);
            b_unique = all_b(1,:);
            
            figure('Name', sprintf('Superfícies C^h(a,b) — Vf=%.2f', vf_target), ...
                   'Position', [50, 50, 1400, 900]);
            
            titles = {'C_{11}', 'C_{22}', 'C_{33} (shear)', ...
                      'C_{12}', 'C_{11}/C_{22}', 'Zener Ratio'};
            
            denom = C11 - C12;
            denom(abs(denom) < 1e-10) = NaN;
            
            data = {C11, C22, C33, C12, C11./C22, 2*C33./denom};
            
            for k = 1:6
                subplot(2,3,k);
                surf(b_unique, a_unique, data{k});
                xlabel('b'); ylabel('a'); 
                title(titles{k}, 'FontSize', 11);
                colorbar; shading interp;
                colormap(jet);
                view(45, 30);
            end

            
            sgtitle(sprintf('C^h(a,b) — Vf=%.2f, holeSize=%.3f', ...
                    vf_target, sqrt(1-vf_target)), 'FontSize', 13);
        end
        function plotMicrosForCurves(obj, a_vals, b_vals, vf_target)
            nA = length(a_vals);
            nB = length(b_vals);
            
            % Linha 1: b=0, a variável
            row1_a = a_vals;
            row1_b = zeros(1, nA);
            
            % Linha 2: a=1, b em 5 pontos
            b_indices = unique([1, round(nB/4), round(nB/2), round(3*nB/4), nB]);
            row2_a = ones(1, length(b_indices));
            row2_b = b_vals(b_indices);
            
            nCols = max(nA, length(b_indices));
            
            figure('Name', sprintf('Micros — Vf=%.2f', vf_target), ...
                   'Position', [100, 50, 200*nCols, 500]);
            
            % Plotar Linha 1
            for idx = 1:nA
                a = row1_a(idx);
                b = row1_b(idx);
                d = (1 + b^2) / a;
                
                obj.latticeVectors = [a, b; b, d];
                obj.defineMesh();
                obj.test = LagrangianFunction.create(obj.baseMesh, 1, 'P1');
                dens = obj.createDensityLevelSet();
                
                subplot(2, nCols, idx);
                rho    = dens.project('P0');
                rhoVal = squeeze(rho.fValues);
                patch('Faces', obj.baseMesh.connec, 'Vertices', obj.baseMesh.coord, ...
                      'FaceVertexCData', rhoVal, 'FaceColor', 'flat', 'EdgeColor', 'none');
                colormap(flipud(gray)); caxis([0 1]);
                axis equal; axis off;
                title(sprintf('a=%.1f, b=0', a), 'FontSize', 8);
            end
            
            % Plotar Linha 2
            for idx = 1:length(b_indices)
                a = row2_a(idx);
                b = row2_b(idx);
                d = (1 + b^2) / a;
                
                obj.latticeVectors = [a, b; b, d];
                obj.defineMesh();
                obj.test = LagrangianFunction.create(obj.baseMesh, 1, 'P1');
                dens = obj.createDensityLevelSet();
                
                subplot(2, nCols, nCols + idx);
                rho    = dens.project('P0');
                rhoVal = squeeze(rho.fValues);
                patch('Faces', obj.baseMesh.connec, 'Vertices', obj.baseMesh.coord, ...
                      'FaceVertexCData', rhoVal, 'FaceColor', 'flat', 'EdgeColor', 'none');
                colormap(flipud(gray)); caxis([0 1]);
                axis equal; axis off;
                title(sprintf('a=1, b=%.1f', b), 'FontSize', 8);
            end
            
            sgtitle(sprintf('Microestruturas — Vf=%.2f\nLinha 1: b=0, a variável | Linha 2: a=1, b variável', vf_target));
        end
        function plotMicrosOneCurve(obj, b_vals, a_vals, vf_target)
            % Seleccionar 6 pontos representativos
            nB      = length(b_vals);
            indices = unique([1, round(nB/5), round(2*nB/5), round(3*nB/5), round(4*nB/5), nB]);
            nPlot   = length(indices);
            
            figure('Name', sprintf('Micros a(b) — Vf=%.2f', vf_target), ...
                   'Position', [100, 50, 200*nPlot, 300]);
            
            for idx = 1:nPlot
                iB = indices(idx);
                b  = b_vals(iB);
                a  = a_vals(iB);
                d  = (1 + b^2) / a;
                
                obj.latticeVectors = [a, b; b, d];
                obj.defineMesh();
                obj.test = LagrangianFunction.create(obj.baseMesh, 1, 'P1');
                dens = obj.createDensityLevelSet();
                
                subplot(1, nPlot, idx);
                rho    = dens.project('P0');
                rhoVal = squeeze(rho.fValues);
                patch('Faces', obj.baseMesh.connec, 'Vertices', obj.baseMesh.coord, ...
                      'FaceVertexCData', rhoVal, 'FaceColor', 'flat', 'EdgeColor', 'none');
                colormap(flipud(gray)); caxis([0 1]);
                axis equal; axis off;
                title(sprintf('b=%.2f\na=%.2f', b, a), 'FontSize', 8);
            end
            
            sgtitle(sprintf('Micros: a(b)=sqrt(1+b²) — Vf=%.2f', vf_target));
        end
        
    end
end