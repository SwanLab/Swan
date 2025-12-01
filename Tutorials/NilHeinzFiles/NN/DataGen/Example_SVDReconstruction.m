% EXAMPLE_SVDRECONSTRUCTION Ejemplo de cómo usar U, S, V para reconstruir T
%
% Este ejemplo demuestra:
% 1. La estructura de U, S, V después del SVD
% 2. Cómo U NO es una función (es una matriz pre-computada)
% 3. Cómo V se predice con NN (fácil de obtener)
% 4. Cómo reconstruir T a partir de U, S, V


fprintf('=== Ejemplo: Estructura de Datos SVD y Reconstrucción de T ===\n\n');

%% ================================================================
%% PASO 1: Generar datos y calcular SVD
%% ================================================================
fprintf('--- PASO 1: Generar datos y calcular SVD ---\n');

radii = 0.1:0.05:0.9;
dataGen = TrainingDataGenerator(radii);
dataGen.generateData(true);  % Generar datos Y calcular SVD

% Obtener resultados SVD
[U, S, V, k] = dataGen.getSVDResults();
fprintf('\nEstructura de datos SVD:\n');
fprintf('  U: [%d × %d] - Modos espaciales (FIJOS, pre-computados)\n', size(U,1), size(U,2));
fprintf('  S: [%d × %d] - Valores singulares (FIJOS, pre-computados)\n', size(S,1), size(S,2));
fprintf('  V: [%d × %d] - Modos paramétricos (varían con r)\n', size(V,1), size(V,2));
fprintf('  k: %d modos retenidos\n', k);

%% ================================================================
%% PASO 2: Explicar que U NO es una función
%% ================================================================
fprintf('\n--- PASO 2: U NO es una función, es una matriz ---\n');

fprintf('\nU es una matriz pre-computada:\n');
fprintf('  - Se calcula UNA VEZ durante el SVD\n');
fprintf('  - Se almacena como matriz [nDofs×8 × k]\n');
fprintf('  - NO necesita inputs (ya está calculado)\n');
fprintf('  - Se usa directamente indexando: U(:,i)\n');

% Ejemplo: Acceder a modos espaciales
fprintf('\nEjemplo de acceso a modos espaciales:\n');
fprintf('  U(:,1) = primer modo espacial [%d × 1]\n', size(U,1));
fprintf('  U(:,2) = segundo modo espacial [%d × 1]\n', size(U,1));
fprintf('  U(:,k) = k-ésimo modo espacial [%d × 1]\n', size(U,1));

% Mostrar valores singulares
sigma = diag(S);
fprintf('\nValores singulares (primeros 5):\n');
for i = 1:min(5, k)
    fprintf('  σ_%d = %.6e\n', i, sigma(i));
end

%% ================================================================
%% PASO 3: Explicar que V se predice con NN (fácil de obtener)
%% ================================================================
fprintf('\n--- PASO 3: V se predice con NN (fácil de obtener) ---\n');

fprintf('\nV se predice con k redes neuronales:\n');
fprintf('  - Entrenar k redes: NN_1, NN_2, ..., NN_k\n');
fprintf('  - Cada red: Input = r (escalar), Output = V(r,i) (escalar)\n');
fprintf('  - Para un radio nuevo r_new:\n');
fprintf('    V_predicted = [NN_1(r_new), NN_2(r_new), ..., NN_k(r_new)]\n');

% Simular predicción de V para un radio nuevo
r_new = 0.35;
fprintf('\nEjemplo: Predicción de V para r_new = %.2f\n', r_new);
fprintf('  (Simulando predicción de NN - en la práctica usarías las redes entrenadas)\n');

% En la práctica, esto sería:
% V_predicted = zeros(k, 1);
% for i = 1:k
%     V_predicted(i) = NN_i.computeOutputValues(r_new);
% end

% Por ahora, interpolamos desde V existente
[~, idx] = min(abs(radii - r_new));
V_predicted = V(idx, :)';  % [k × 1]
fprintf('  V_predicted = [%.6f, %.6f, ..., %.6f]'' [%d × 1]\n', ...
    V_predicted(1), V_predicted(2), V_predicted(end), length(V_predicted));

%% ================================================================
%% PASO 4: Reconstruir T a partir de U, S, V
%% ================================================================
fprintf('\n--- PASO 4: Reconstruir T a partir de U, S, V ---\n');

% Reconstrucción: T = U × S × V'
fprintf('\nFórmula de reconstrucción:\n');
fprintf('  T_vec = U × S × V_predicted\n');
fprintf('  T = reshape(T_vec, [nDofs, 8])\n');

% Obtener nDofs del primer resultado
firstResult = dataGen.allResults{1};
nDofs = size(firstResult.u, 1);

% Reconstruir T
T_vec = U * S * V_predicted;  % [nDofs×8 × 1]
T_reconstructed = reshape(T_vec, [nDofs, 8]);  % [nDofs × 8]

fprintf('\nResultado:\n');
fprintf('  T_vec: [%d × 1] (vectorizado)\n', size(T_vec,1));
fprintf('  T_reconstructed: [%d × 8] (forma original)\n', size(T_reconstructed,1), size(T_reconstructed,2));

% Comparar con T original
T_original = firstResult.u;  % Para el radio más cercano
reconstruction_error = norm(T_reconstructed - T_original, 'fro') / norm(T_original, 'fro');
fprintf('\nError de reconstrucción (comparado con T original): %.2e\n', reconstruction_error);

%% ================================================================
%% PASO 5: Visualización de la estructura
%% ================================================================
fprintf('\n--- PASO 5: Visualización de la estructura ---\n');

figure('Position', [100, 100, 1400, 800]);

% Subplot 1: Valores singulares
subplot(2, 3, 1);
semilogy(1:k, sigma(1:k), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Modo i');
ylabel('Valor singular σ_i');
title('Valores Singulares');
grid on;

% Subplot 2: Primeros modos paramétricos V
subplot(2, 3, 2);
nModesToPlot = min(5, k);
hold on;
for i = 1:nModesToPlot
    plot(radii, V(:, i), '-', 'LineWidth', 1.5, 'DisplayName', sprintf('V_%d', i));
end
xlabel('Radio r');
ylabel('Coeficiente V_i(r)');
title('Modos Paramétricos V (primeros 5)');
legend('Location', 'best');
grid on;

% Subplot 3: Predicción de V para r_new
subplot(2, 3, 3);
bar(1:k, V_predicted, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Modo i');
ylabel('V_i(r_{new})');
title(sprintf('V Predicho para r = %.2f', r_new));
grid on;

% Subplot 4: Estructura de U (primeros valores de cada modo)
subplot(2, 3, 4);
nSamples = min(100, size(U,1));
imagesc(U(1:nSamples, 1:min(5,k)));
colorbar;
xlabel('Modo i');
ylabel('Componente (primeras 100)');
title('Estructura de U (primeros 100 componentes)');
colormap('jet');

% Subplot 5: Comparación T original vs reconstruido (primer modo)
subplot(2, 3, 5);
plot(1:nDofs, T_original(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'T original');
hold on;
plot(1:nDofs, T_reconstructed(:,1), 'r--', 'LineWidth', 1.5, 'DisplayName', 'T reconstruido');
xlabel('Grado de libertad');
ylabel('Desplazamiento');
title('Comparación T (modo 1)');
legend('Location', 'best');
grid on;

% Subplot 6: Error de reconstrucción por modo
subplot(2, 3, 6);
errors_per_mode = zeros(8, 1);
for mode = 1:8
    errors_per_mode(mode) = norm(T_reconstructed(:,mode) - T_original(:,mode)) / ...
                            norm(T_original(:,mode));
end
bar(1:8, errors_per_mode, 'FaceColor', [0.8 0.2 0.2]);
xlabel('Modo grueso');
ylabel('Error relativo');
title('Error de Reconstrucción por Modo');
grid on;

sgtitle('Estructura de Datos SVD y Reconstrucción de T', 'FontSize', 14, 'FontWeight', 'bold');

%% ================================================================
%% RESUMEN
%% ================================================================
fprintf('\n=== RESUMEN ===\n');
fprintf('\n1. U (Modos Espaciales):\n');
fprintf('   - Tipo: Matriz pre-computada [%d × %d]\n', size(U,1), size(U,2));
fprintf('   - NO es una función, es una matriz fija\n');
fprintf('   - Input: Índice del modo i (1 a k)\n');
fprintf('   - Output: Vector U(:,i) [%d × 1]\n', size(U,1));
fprintf('   - NO depende de r (fijo)\n');

fprintf('\n2. S (Valores Singulares):\n');
fprintf('   - Tipo: Matriz diagonal [%d × %d]\n', size(S,1), size(S,2));
fprintf('   - Fijo, pre-computado\n');
fprintf('   - NO depende de r (fijo)\n');

fprintf('\n3. V (Modos Paramétricos):\n');
fprintf('   - Tipo: Se predice con redes neuronales\n');
fprintf('   - Input: Radio r (escalar)\n');
fprintf('   - Output: Vector V(r) [%d × 1]\n', k);
fprintf('   - SÍ depende de r (se predice con NN)\n');
fprintf('   - FÁCIL de obtener: solo necesitas k redes neuronales\n');

fprintf('\n4. T (Reconstruido):\n');
fprintf('   - Tipo: Matriz reconstruida [%d × 8]\n', nDofs);
fprintf('   - Fórmula: T = reshape(U × S × V(r), [%d, 8])\n', nDofs);
fprintf('   - Se reconstruye cuando se necesita\n');

fprintf('\n✓ La estructura permite:\n');
fprintf('   - U y S se calculan UNA VEZ (costoso)\n');
fprintf('   - V se predice RÁPIDO con NN (barato)\n');
fprintf('   - T se reconstruye RÁPIDO (multiplicación de matrices)\n');

