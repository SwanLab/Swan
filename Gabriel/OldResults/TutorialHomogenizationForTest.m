% %% TEST_RHOMBOID_ALPHA - usando TutorialHomogenization
% 
% clear; close all; clc;
% 
% c1 = 1.0;
% Area = 1.0;
% alpha_vals = 30:5:150;
% 
% nAlpha = length(alpha_vals);
% 
% C11_vals = zeros(nAlpha,1);
% C22_vals = zeros(nAlpha,1);
% C12_vals = zeros(nAlpha,1);
% C33_vals = zeros(nAlpha,1);
% vf_vals  = zeros(nAlpha,1);
% 
% for i = 1:nAlpha
% 
%     alpha = alpha_vals(i);
%     c2 = Area / (c1*sind(alpha));
% 
%     fprintf('α = %4d° → c2 = %.4f ', alpha, c2);
% 
%     try
%         TH = TutorialHomogenization(c1,c2,alpha,40);
%         [C,vf] = TH.runAtVolume(0.5);
% 
%         C11_vals(i) = C(1,1,1,1);
%         C22_vals(i) = C(2,2,2,2);
%         C12_vals(i) = C(1,1,2,2);
%         C33_vals(i) = C(1,2,1,2);
%         vf_vals(i)  = vf;
% 
%         fprintf('→ vf=%.3f\n', vf);
% 
%     catch ME
%         fprintf('→ ERROR: %s\n', ME.message);
%         C11_vals(i) = NaN;
%         C22_vals(i) = NaN;
%         C12_vals(i) = NaN;
%         C33_vals(i) = NaN;
%         vf_vals(i)  = NaN;
%     end
% end
% 
% % 🔥 EXATAMENTE o mesmo plot que você já tinha
% plotAlphaResults(alpha_vals, ...
%                  C11_vals, C22_vals, ...
%                  C12_vals, C33_vals, ...
%                  vf_vals);
% %% =====================================================
% %% FUNÇÃO DE PLOT
% %% =====================================================
% 
% function plotAlphaResults(alpha_vals, C11, C22, C12, C33, vf)
% 
%     colors = [
%         0.2667, 0.4667, 0.6667;  % blue
%         0.9333, 0.4000, 0.4667;  % red
%         0.4000, 0.7333, 0.4000;  % green
%         0.9333, 0.6667, 0.2000;  % orange
%         0.6000, 0.4000, 0.6667;  % purple
%     ];
% 
%     figure('Color','white','Position',[100 100 1400 1000])
% 
%     subplot(2,2,1)
%     plot(alpha_vals,C11,'o-','LineWidth',1.5,'Color',colors(1,:))
%     xlabel('\alpha (deg)')
%     ylabel('C_{11}')
%     grid on
% 
%     subplot(2,2,2)
%     plot(alpha_vals,C22,'s-','LineWidth',1.5,'Color',colors(2,:))
%     xlabel('\alpha (deg)')
%     ylabel('C_{22}')
%     grid on
% 
%     subplot(2,2,3)
%     plot(alpha_vals,C33,'d-','LineWidth',1.5,'Color',colors(3,:))
%     xlabel('\alpha (deg)')
%     ylabel('C_{33}')
%     grid on
% 
%     subplot(2,2,4)
%     plot(alpha_vals,C12,'^-','LineWidth',1.5,'Color',colors(4,:))
%     xlabel('\alpha (deg)')
%     ylabel('C_{12}')
%     grid on
% 
%     sgtitle('Homogenized Tensor vs Angle \alpha (vf ≈ 0.5)')
% 
%     % Anisotropy ratio
%     figure('Color','white')
%     plot(alpha_vals,C11./C22,'*-','LineWidth',2,'Color',colors(5,:))
%     xlabel('\alpha (deg)')
%     ylabel('C_{11}/C_{22}')
%     grid on
%     yline(1,'--','Isotropy')
% 
%     % Volume fraction check
%     figure('Color','white')
%     plot(alpha_vals,vf,'o-','LineWidth',1.5)
%     xlabel('\alpha (deg)')
%     ylabel('vf')
%     grid on
%     yline(0.5,'--','Target')
% end

%% TEST_RHOMBOID_ALPHA - usando TutorialHomogenization

clear; close all; clc;

c1 = 1.0;
Area = 1.0;
b_vals = linspace(0,1,40);

nAlpha = length(b_vals);

C11_vals = zeros(nAlpha,1);
C22_vals = zeros(nAlpha,1);
C12_vals = zeros(nAlpha,1);
C33_vals = zeros(nAlpha,1);
vf_vals  = zeros(nAlpha,1);

C1_vals = zeros(nAlpha,1);
C2_vals = zeros(nAlpha,1);
C3_vals = zeros(nAlpha,1);
C4_vals = zeros(nAlpha,1);
C5_vals = zeros(nAlpha,1);

A_vals = zeros(nAlpha,1);
C11_norm = zeros(nAlpha,1);
C22_norm = zeros(nAlpha,1);

for i = 1:nAlpha
    
    b = b_vals(i);
    
    
    % evitar singularidade
    if abs(b) < 1e-6
        theta = 90;   % caso limite
    else
        theta = atan(1/b) * 180/pi;
    end
    
    c2 = abs(Area / (c1*sind(theta)));
    
    fprintf('b = %6.3f → θ = %6.2f° → c2 = %.4f ', b, theta, c2);
    
    try
        TH = TutorialHomogenization(c1,c2,theta,40);
        [C,vf] = TH.runAtVolume(0.5);
        C_all(:,:,:,:,i) = C;
        
        C11_vals(i) = C(1,1,1,1);
        C22_vals(i) = C(2,2,2,2);
        C12_vals(i) = C(1,1,2,2);
        C33_vals(i) = C(1,2,1,2);
        vf_vals(i)  = vf;

        C11 = C11_vals(i);
        C22 = C22_vals(i);
        C12 = C12_vals(i);
        C33 = C33_vals(i);
        
        C1_vals(i) = (1/8)*(3*(C11+C22) + 2*(C12 + 2*C33));
        C2_vals(i) = (1/2)*(C11 - C22);
        C3_vals(i) = (1/8)*((C11+C22) - 2*(C12 + 2*C33));
        C4_vals(i) = (1/8)*((C11+C22) + 2*(3*C12 - 2*C33));
        C5_vals(i) = (1/8)*((C11+C22) - 2*(C12 - 2*C33));

        A_vals(i) = 2*C33_vals(i) / (C11_vals(i) - C12_vals(i));
        
        fprintf('→ vf=%.3f\n', vf);
        
    catch ME
        fprintf('→ ERROR: %s\n', ME.message);
        
        C11_vals(i) = NaN;
        C22_vals(i) = NaN;
        C12_vals(i) = NaN;
        C33_vals(i) = NaN;
        vf_vals(i)  = NaN;


    end
end
% C11_norm = C11_vals / C11_vals(1);
% C22_norm = C22_vals / C22_vals(1);
ref = find(~isnan(C11_vals),1);

C11_norm = C11_vals / C11_vals(ref);
C22_norm = C22_vals / C22_vals(ref);



% 🔥 EXATAMENTE o mesmo plot que você já tinha
plotAlphaResults(b_vals, ...
                 C11_vals, C22_vals, ...
                 C12_vals, C33_vals, ...
                 vf_vals,A_vals, C11_norm, C22_norm,C1_vals,C2_vals,C3_vals,C4_vals,C5_vals);
%% =====================================================
%% FUNÇÃO DE PLOT
%% =====================================================

function plotAlphaResults(b_vals, C11, C22, C12, C33, vf, A_vals, C11_norm, C22_norm,C1_vals,C2_vals,C3_vals,C4_vals,C5_vals)

    colors = [
        0.2667, 0.4667, 0.6667;  % blue
        0.9333, 0.4000, 0.4667;  % red
        0.4000, 0.7333, 0.4000;  % green
        0.9333, 0.6667, 0.2000;  % orange
        0.6000, 0.4000, 0.6667;  % purple
    ];

    figure('Color','white','Position',[100 100 1400 1000])

    subplot(2,2,1)
    plot(b_vals,C11,'o-','LineWidth',1.5,'Color',colors(1,:))
    xlabel('\alpha (deg)')
    ylabel('C_{11}')
    grid on

    subplot(2,2,2)
    plot(b_vals,C22,'s-','LineWidth',1.5,'Color',colors(2,:))
    xlabel('\alpha (deg)')
    ylabel('C_{22}')
    grid on

    subplot(2,2,3)
    plot(b_vals,C33,'d-','LineWidth',1.5,'Color',colors(3,:))
    xlabel('\alpha (deg)')
    ylabel('C_{33}')
    grid on

    subplot(2,2,4)
    plot(b_vals,C12,'^-','LineWidth',1.5,'Color',colors(4,:))
    xlabel('\alpha (deg)')
    ylabel('C_{12}')
    grid on

    sgtitle('Homogenized Tensor vs Angle \alpha (vf ≈ 0.5)')

    % Anisotropy ratio
    figure('Color','white')
    plot(b_vals,C11./C22,'*-','LineWidth',2,'Color',colors(5,:))
    xlabel('\alpha (deg)')
    ylabel('C_{11}/C_{22}')
    grid on
    yline(1,'--','Isotropy')

    % Volume fraction check
    figure('Color','white')
    plot(b_vals,vf,'o-','LineWidth',1.5)
    xlabel('\alpha (deg)')
    ylabel('vf')
    grid on
    yline(0.5,'--','Target')

    figure('Color','white')
    plot(b_vals, A_vals,'o-','LineWidth',2)
    xlabel('b')
    ylabel('Zener ratio A(b)')
    grid on
    yline(1,'--','Isotropy')
    title('Zener anisotropy ratio')
    figure('Color','white')
    
    plot(b_vals,C11_norm,'o-','LineWidth',1.5); hold on
    plot(b_vals,C22_norm,'s-','LineWidth',1.5)
    
    xlabel('b')
    ylabel('Normalized stiffness')
    legend('C11/C11(0)','C22/C22(0)')
    grid on
    
    title('Normalized elastic components')

    figure('Color','white')

    subplot(3,2,1)
    plot(b_vals,C1_vals,'LineWidth',1.5)
    title('C1')
    grid on
    
    subplot(3,2,2)
    plot(b_vals,C2_vals,'LineWidth',1.5)
    title('C2 (anisotropy)')
    grid on
    
    subplot(3,2,3)
    plot(b_vals,C3_vals,'LineWidth',1.5)
    title('C3')
    grid on
    
    subplot(3,2,4)
    plot(b_vals,C4_vals,'LineWidth',1.5)
    title('C4')
    grid on
    
    subplot(3,2,5)
    plot(b_vals,C5_vals,'LineWidth',1.5)
    title('C5')
    grid on
    
    sgtitle('Elastic invariants vs b')


end

figure
pax = polaraxes;
hold(pax,'on')

idx_list = [1, round(nAlpha/2), nAlpha];

for k = 1:length(idx_list)
    
    i = idx_list(k);
    
    Ctest = C_all(:,:,:,:,i);
    
    [theta, Etheta] = computePolarE(Ctest);
    
    polarplot(pax, theta, Etheta,'LineWidth',2)
end

legend('b min','b mid','b max')
title('Directional Young modulus E(\theta)')

function [theta, Etheta] = computePolarE(C)

    % --- tensor para Voigt ---
    Cvoigt = zeros(3,3);
    Cvoigt(1,1) = C(1,1,1,1);
    Cvoigt(1,2) = C(1,1,2,2);
    Cvoigt(2,1) = C(2,2,1,1);
    Cvoigt(2,2) = C(2,2,2,2);
    Cvoigt(3,3) = C(1,2,1,2);

    % --- compliance ---
    S = pinv(Cvoigt);

    % --- ângulos ---
    theta = linspace(0,2*pi,200);
    Etheta = zeros(size(theta));

    for i = 1:length(theta)
        c = cos(theta(i));
        s = sin(theta(i));

        n = [c^2;
             s^2;
             2*c*s];

        Etheta(i) = 1 / (n' * S * n);
    end
end
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
        pnorm
        monitoring
        Mmass
        fixedHoleSize
        currentB
        latticeVectors
        lastFem
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
            obj.computeHoleParams();
            obj.compute();
            obj.plot();
            obj.printTensorAtVolume(0.5);
            obj.verifyTensor();
        end 
        
        %% ========== MÉTODOS PÚBLICOS DE VISUALIZAÇÃO ==========

       function plotDeformedStates(obj, b_idx, scale)
    % Plota as deformadas para os três casos de carga
    % b_idx: índice do valor de b (opcional, padrão = último)
    % scale: fator de escala para deformação (opcional, padrão = 1)
    %        scale = 1 → deformação real (pode ser muito pequena)
    %        scale = 10 ou 20 → amplificada para visualização
    
            if nargin < 2 || isempty(b_idx)
                b_idx = length(obj.paramHole);
            end
            if nargin < 3 || isempty(scale)
                scale = 1;
            end
                    
            b_val = double(obj.paramHole(b_idx));
            
            % Recalcular o problema para este b específico
            v1 = [1, 0];
            v2 = [b_val, 1];
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
            
            figure('Name', sprintf('Deformed States - b = %.4f', b_val), 'Position', [100, 100, 1400, 500]);
            
            for iBasis = 1:3
                subplot(1,3,iBasis);
                
                uFlucFun = fem.uFluc{iBasis};
                uFluc = uFlucFun.fValues;
                
                switch iBasis
                    case 1
                        E_macro = [1, 0; 0, 0];
                        title_str = 'Tensão Horizontal (ε_{xx}=1)';
                    case 2
                        E_macro = [0, 0; 0, 1];
                        title_str = 'Tensão Vertical (ε_{yy}=1)';
                    case 3
                        E_macro = [0, 0.5; 0.5, 0];
                        title_str = 'Cisalhamento (ε_{xy}=1)';
                end
                
                coord = obj.baseMesh.coord;
                u_macro = (E_macro * coord')';
                u_total = u_macro + uFluc;
                coord_def = coord + scale * u_total;
                
                trisurf(obj.baseMesh.connec, coord_def(:,1), coord_def(:,2), zeros(size(coord_def,1),1), ...
                        'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5);
                view(2);
                axis equal;
                title(title_str);
                xlabel('x'); ylabel('y');
                
                hold on;
                % Destacar vértices
                plot(coord_def(1:4,1), coord_def(1:4,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                
                % Destacar nós master e slave
                masterNodes = obj.masterSlave(:,1);
                slaveNodes = obj.masterSlave(:,2);
                plot(coord_def(masterNodes,1), coord_def(masterNodes,2), 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
                plot(coord_def(slaveNodes,1), coord_def(slaveNodes,2), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
                hold off;
            end
            
            if scale == 1
                sgtitle(sprintf('Deformadas para b = %.4f (escala real = 1)', b_val));
            else
                sgtitle(sprintf('Deformadas para b = %.4f (escala = %d)', b_val, scale));
            end
        end

        function checkPeriodicity(obj, b_idx)
            % Verifica a condição periódica: deslocamentos flutuantes devem ser iguais nos pares master-slave
            if nargin < 2 || isempty(b_idx)
                b_idx = length(obj.paramHole);
            end
            
            b_val = double(obj.paramHole(b_idx)); 

            v1 = [1, 0];
            v2 = [b_val, 1];
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
            
            figure('Name', sprintf('Periodicity Check - b = %.4f', b_val), 'Position', [100, 100, 1200, 400]);
            
            masterNodes = obj.masterSlave(:,1);
            slaveNodes = obj.masterSlave(:,2);
            
            % Ordenar por posição y para visualização mais limpa
            [~, sortIdx] = sort(obj.baseMesh.coord(masterNodes, 2));
            
            for iBasis = 1:3
                subplot(1,3,iBasis);
                
                uFlucFun = fem.uFluc{iBasis};
                uFluc = uFlucFun.fValues;
                
                % Deslocamentos flutuantes nos nós master e slave
                uFluc_master = uFluc(masterNodes, :);
                uFluc_slave  = uFluc(slaveNodes, :);
                
                % A diferença DEVE ser zero (ou muito próxima) para condições periódicas
                diff_flutu = uFluc_slave - uFluc_master;
                err = vecnorm(diff_flutu, 2, 2);
                
                switch iBasis
                    case 1
                        titulo = 'Base 1: Tração Horizontal';
                    case 2
                        titulo = 'Base 2: Tração Vertical';
                    case 3
                        titulo = 'Base 3: Cisalhamento';
                end
                
                % Plotar erro para cada par
                plot(1:length(err), err(sortIdx), 'b.', 'MarkerSize', 6);
                hold on
                yline(1e-10, 'r--', 'Tolerância 1e-10', 'LineWidth', 1.5);
                yline(1e-12, 'g--', 'Tolerância 1e-12', 'LineWidth', 1);
                
                xlabel('Par master-slave'); 
                ylabel('||u^{fluc}_{slave} - u^{fluc}_{master}||');
                title(sprintf('%s\nErro máximo: %.2e', titulo, max(err)));
                set(gca, 'YScale', 'log');
                grid on;
                legend('Erro por par', 'Tol 1e-10', 'Tol 1e-12', 'Location', 'best');
                xlim([1, length(err)]);
                hold off
            end
            
            sgtitle(sprintf('Verificação de Periodicidade - b = %.4f\n(Deslocamentos flutuantes devem ser iguais nos pares master-slave)', b_val));
        end

        function verifyTensor(obj, idx)
            if nargin < 2
                idx = length(obj.paramHole);
            end
            
            b_val = obj.paramHole(idx);
            C = obj.Chomog(:,:,:,:,idx);
            vf = obj.volFrac(idx);
            
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
            
            fprintf('\n══════════════════════════════════════════════════════════════\n');
            fprintf('              VERIFICAÇÃO DO TENSOR HOMOGENEIZADO              \n');
            fprintf('══════════════════════════════════════════════════════════════\n');
            fprintf('b = %.6f\n', b_val);
            fprintf('Volume fraction = %.4f\n', vf);
            fprintf('──────────────────────────────────────────────────────────────\n');
            
            % Simetria
            fprintf('\n1. SIMETRIA DO TENSOR:\n');
            sym_12 = abs(Cvoigt(1,2) - Cvoigt(2,1));
            sym_13 = abs(Cvoigt(1,3) - Cvoigt(3,1));
            sym_23 = abs(Cvoigt(2,3) - Cvoigt(3,2));
            
            fprintf('   C12 = %.6f, C21 = %.6f → dif = %.2e %s\n', ...
                    Cvoigt(1,2), Cvoigt(2,1), sym_12, obj.check(sym_12 < 1e-6));
            fprintf('   C13 = %.6f, C31 = %.6f → dif = %.2e %s\n', ...
                    Cvoigt(1,3), Cvoigt(3,1), sym_13, obj.check(sym_13 < 1e-6));
            fprintf('   C23 = %.6f, C32 = %.6f → dif = %.2e %s\n', ...
                    Cvoigt(2,3), Cvoigt(3,2), sym_23, obj.check(sym_23 < 1e-6));
            
            % Estabilidade
            fprintf('\n2. ESTABILIDADE (positividade definida):\n');
            det1 = Cvoigt(1,1);
            det2 = Cvoigt(1,1)*Cvoigt(2,2) - Cvoigt(1,2)^2;
            det3 = det(Cvoigt);
            
            fprintf('   C11 = %.6f > 0 → %s\n', det1, obj.check(det1 > 0));
            fprintf('   C11*C22 - C12^2 = %.6f > 0 → %s\n', det2, obj.check(det2 > 0));
            fprintf('   det(C) = %.6f > 0 → %s\n', det3, obj.check(det3 > 0));
            
            eigvals = eig(Cvoigt);
            fprintf('\n   Autovalores:\n');
            fprintf('   λ1 = %.6f\n', eigvals(1));
            fprintf('   λ2 = %.6f\n', eigvals(2));
            fprintf('   λ3 = %.6f → %s\n', eigvals(3), obj.check(all(eigvals > 0)));
            
            % Consistência energética
            fprintf('\n3. CONSISTÊNCIA ENERGÉTICA (ε:C:ε > 0):\n');
            eps_x = [1; 0; 0];
            E_x = eps_x' * Cvoigt * eps_x;
            fprintf('   Tração uniaxial em x: E = %.6f → %s\n', E_x, obj.check(E_x > 0));
            
            eps_y = [0; 1; 0];
            E_y = eps_y' * Cvoigt * eps_y;
            fprintf('   Tração uniaxial em y: E = %.6f → %s\n', E_y, obj.check(E_y > 0));
            
            eps_s = [0; 0; 1];
            E_s = eps_s' * Cvoigt * eps_s;
            fprintf('   Cisalhamento puro: E = %.6f → %s\n', E_s, obj.check(E_s > 0));
            
            % Módulos efetivos
            fprintf('\n4. MÓDULOS EFETIVOS:\n');
            K = (Cvoigt(1,1) + Cvoigt(2,2) + 2*Cvoigt(1,2))/4;
            G = Cvoigt(3,3);
            Eeff = G*(3*K + G)/(K + G);
            nueff = (K - G)/(K + G);
            
            fprintf('   Bulk modulus efetivo K = %.6f\n', K);
            fprintf('   Shear modulus efetivo G = %.6f\n', G);
            fprintf('   Young modulus efetivo E = %.6f\n', Eeff);
            fprintf('   Poisson ratio efetivo ν = %.6f\n', nueff);
            
            % Anisotropia
            fprintf('\n5. ANISOTROPIA:\n');
            A = 2*Cvoigt(3,3) / (Cvoigt(1,1) - Cvoigt(1,2));
            fprintf('   Zener ratio A = %.6f\n', A);
            fprintf('   C11/C22 = %.6f\n', Cvoigt(1,1)/Cvoigt(2,2));
            
            fprintf('\n──────────────────────────────────────────────────────────────\n');
            
            all_checks = (sym_12 < 1e-6) && (sym_13 < 1e-6) && (sym_23 < 1e-6) && ...
                         (det1 > 0) && (det2 > 0) && (det3 > 0) && ...
                         (E_x > 0) && (E_y > 0) && (E_s > 0);
            
            if all_checks
                fprintf('✅ RESULTADO: TENSOR FISICAMENTE CONSISTENTE\n');
            else
                fprintf('❌ RESULTADO: TENSOR COM INCONSISTÊNCIAS\n');
            end
            fprintf('══════════════════════════════════════════════════════════════\n\n');
        end

        function plot(obj)
            tiledlayout(1,3)
            nexttile
            hold on
            plot(obj.paramHole,squeeze(obj.Chomog(1,1,1,1,:)),'LineStyle','none','Marker','o')
            nexttile
            hold on
            plot(obj.paramHole,squeeze(obj.Chomog(2,2,2,2,:)),'LineStyle','none','Marker','o')
            nexttile
            hold on
            plot(obj.paramHole,squeeze(obj.Chomog(1,2,1,2,:)),'LineStyle','none','Marker','o')
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.E          = 1;
            obj.nu         = 0.3;
            obj.meshType   = 'Square';
            obj.meshN      = 80;
            obj.holeType   = 'Square';
            obj.pnorm      = 'Inf';
            obj.nSteps     = 80;
            obj.monitoring = true;
            obj.fixedHoleSize = 1e-10;
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
            obj.paramHole = linspace(1e-6, obj.maxParam, obj.nSteps);
        end

        function compute(obj)
            nComb = length(obj.paramHole);
            mat = zeros(2,2,2,2,nComb);
            volF = zeros(1,nComb);
            
            fprintf('========================================\n');
            fprintf('b varia de %.2e a %.2f\n', obj.paramHole(1), obj.paramHole(end));
            fprintf('Tamanho do furo FIXO: %.2f\n', obj.fixedHoleSize);
            fprintf('========================================\n\n');
            
            for i = 1:nComb
                b_val = obj.paramHole(i);
                obj.currentB = b_val;             
                
                v1 = [1, 0];
                v2 = [b_val, 1];
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
            fem.plotMicroFields(2,0);
            
            obj.lastFem = fem;

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
            fprintf('Volume encontrado:  %.4f (b = %.4f)\n', obj.volFrac(idx), obj.paramHole(idx));
        
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
            title(['Tensor homogenizado - b = ', num2str(obj.paramHole(idx))])
        end
        
    end
    
    methods (Access = private, Static)
        
        function str = check(condition, true_msg, false_msg)
            if nargin < 2
                true_msg = '✓ OK';
            end
            if nargin < 3
                false_msg = '✗ PROBLEMA';
            end
            
            if condition
                str = true_msg;
            else
                str = false_msg;
            end
        end
        
    end
    
end