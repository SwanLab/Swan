# Análisis Descomposición en Valores Singulares para Matrices T

## 1. Objetivo del Código

Se aplica **Singular Value Decomposition (SVD)** a múltiples matrices `T` (una por cada radio) para:
- **Identificar modos dominantes** en la variación de `T` respecto al radio
- **Reducir la dimensionalidad** del problema
- **Comprender la estructura** de cómo `T` depende del radio `r`

## 2. Estructura de Datos

### 2.1. Matriz T Individual

Cada archivo `.mat` contiene:
- **`T`**: Matriz `[nDofs × 8]` donde:
  - `nDofs = mesh.nnodes × mesh.ndim = mesh.nnodes × 2` (grados de libertad)
  - 8 columnas = 8 modos gruesos
  - Cada columna `T(:,j)` = desplazamientos finos para el modo grueso `j`
- **`mesh`**: Información de la malla (coordenadas, conectividad)

**Ejemplo**: Para una malla 50×50 elementos:
- `nnodes ≈ 2601` nodos
- `nDofs = 2601 × 2 = 5202`
- `T` tiene tamaño `[5202 × 8]`

### 2.2. Construcción de T_SVD

```matlab
r = 0:0.02:0.999;  % 50 radios (0, 0.02, 0.04, ..., 0.98)
```

**Proceso**:
1. Para cada radio `r(i)`:
   - Carga `T` y `mesh` del archivo correspondiente
   - Vectoriza `T` → `T(:)` → vector columna de tamaño `[nDofs × 8, 1]`
   - Almacena en columna `i` de `T_SVD`

2. **Estructura de T_SVD**:
   ```
   T_SVD = [T(:)_r1, T(:)_r2, T(:)_r3, ..., T(:)_r50]
   ```
   - **Filas**: `nDofs × 8 = 5202 × 8 = 41,616` (todas las componentes de T)
   - **Columnas**: `size(r,2) = 50` (una por cada radio)
   - **Dimensiones**: `[41,616 × 50]`

**Interpretación física**:
- Cada **columna** de `T_SVD` = vectorización completa de `T` para un radio específico
- Cada **fila** de `T_SVD` = cómo varía una componente específica de `T` al cambiar el radio

## 3. Descomposición SVD

```matlab
[U,S,V] = svd(T_SVD,'econ');
```

### 3.1. Descomposición Matemática

**SVD descompone**:
```
T_SVD = U × S × V'
```

Donde:
- **`U`**: `[41,616 × min(41,616, 50)] = [41,616 × 50]`
  - Columnas = **modos espaciales** (patrones de desplazamiento en el espacio)
  - Cada columna `U(:,i)` = un modo espacial dominante
  - Ordenados por importancia (primeros modos = más importantes)
  
- **`S`**: `[50 × 50]` (diagonal)
  - Valores singulares: `σ₁ ≥ σ₂ ≥ ... ≥ σ₅₀ ≥ 0`
  - `σᵢ` = importancia del modo `i`
  - Decaimiento rápido → pocos modos dominantes
  
- **`V`**: `[50 × 50]`
  - Columnas = **modos paramétricos** (dependencia respecto al radio)
  - Cada columna `V(:,i)` = cómo varía el modo `i` con el radio `r`
  - `V(j,i)` = coeficiente del modo `i` para el radio `r(j)`

### 3.2. Interpretación Física

**Reconstrucción de T para un radio r(j)**:
```
T_SVD(:,j) ≈ Σᵢ₌₁ᵏ σᵢ × U(:,i) × V(j,i)
```

Donde `k` es el número de modos retenidos (típicamente `k << 50`).

**Significado**:
- **U(:,i)**: Patrón espacial del modo `i` (campo de desplazamientos)
- **V(j,i)**: Peso del modo `i` para el radio `r(j)`
- **σᵢ**: Importancia global del modo `i`

## 4. Visualización de Resultados

### 4.1. Visualización de V (Modos Paramétricos)

```matlab
plot(r, V(:,idx), 'LineWidth', 1.5);
```

**Qué muestra**:
- **Eje X**: Radio `r` (0 a 0.98)
- **Eje Y**: Coeficiente `V(j,idx)` del modo `idx` para cada radio
- **Interpretación**: Cómo varía la importancia del modo `idx` con el radio

**Ejemplo**:
- Si `V(:,1)` es constante → el modo 1 es igualmente importante para todos los radios
- Si `V(:,2)` tiene un pico en `r=0.5` → el modo 2 es más importante para radios intermedios

### 4.2. Visualización de S (Valores Singulares)

```matlab
plot(log(diag(S)), 'LineWidth', 1.5);
```

**Qué muestra**:
- **Eje X**: Índice del modo (1, 2, 3, ..., 50)
- **Eje Y**: `log(σᵢ)` (logaritmo del valor singular)
- **Interpretación**: 
  - Decaimiento rápido → pocos modos dominantes
  - Si `log(σ₁₀) << log(σ₁)` → los primeros 10 modos capturan la mayoría de la información

**Criterio de truncamiento**:
- Retener modos hasta que `σᵢ/σ₁ > tol` (ej: `tol = 1e-6`)
- Típicamente: `k ≈ 10-20` modos son suficientes

### 4.3. Visualización de U (Modos Espaciales)

```matlab
plot(U(:,i));
```

**Qué muestra**:
- **Eje X**: Índice del grado de libertad (1 a 41,616)
- **Eje Y**: Valor del modo espacial `U(:,i)` en ese grado de libertad
- **Interpretación**: Patrón espacial del modo `i` (campo de desplazamientos)

**Nota**: El código comentado muestra cómo visualizar `U` como función sobre la malla:
```matlab
s.mesh = mesh;
s.order = 'P1';
s.fValues = U(:,i);
Ufun = LagrangianFunction(s);
Ufun.plot();  % Visualiza el modo como campo sobre la malla
```

## 5. Aplicación: Reducción de Dimensionalidad

### 5.1. Reconstrucción Aproximada

**Reconstrucción completa** (50 modos):
```
T_SVD(:,j) = U × S × V(j,:)'
```

**Reconstrucción truncada** (k modos, k < 50):
```
T_SVD_approx(:,j) = U(:,1:k) × S(1:k,1:k) × V(j,1:k)'
```

**Error de aproximación**:
```
Error = ||T_SVD(:,j) - T_SVD_approx(:,j)|| / ||T_SVD(:,j)||
```

Si `σₖ₊₁/σ₁ << 1`, el error es pequeño.

### 5.2. Ventajas para Entrenamiento de Redes Neuronales

**Opción 3 sin SVD**:
- Entrenar 8 redes neuronales
- Cada red: Input `[r, x, y]` → Output `[Tx_i, Ty_i]`
- Total: `8 × (nDofs × 2)` = `8 × 10,404` = `83,232` parámetros de salida

**Opción 3 con SVD**:
- Entrenar `k` redes neuronales (k ≈ 10-20)
- Cada red: Input `[r]` → Output `V(j,i)` (coeficiente del modo `i` para radio `r(j)`)
- Reconstrucción: `T = U(:,1:k) × S(1:k,1:k) × V_predicted'`
- Total: `k` parámetros de salida (mucho menor)

**Reducción**:
- De `83,232` parámetros → `k ≈ 10-20` parámetros
- Reducción de ~4000× a ~8000×

## 6. Flujo Completo del Proceso

```
1. Generación de datos 
   └─> Para cada r: T(r) [nDofs × 8]
   
2. SVD 
   └─> T_SVD = [T(:)_r1, T(:)_r2, ..., T(:)_r50]
   └─> [U, S, V] = svd(T_SVD)
   └─> U: modos espaciales [nDofs×8 × k]
   └─> S: valores singulares [k × k]
   └─> V: modos paramétricos [50 × k]
   
3. Entrenamiento de Red Neuronal (Option 3 con SVD)
   └─> Entrenar k redes: r → V(r,i) para i=1..k
   
4. Predicción
   └─> Dado r nuevo:
       ├─> Predecir V_predicted = [V(r,1), V(r,2), ..., V(r,k)]
       └─> Reconstruir: T_predicted = U(:,1:k) × S(1:k,1:k) × V_predicted'
```

## 7. Comparación: Con vs Sin SVD

| Aspecto | Sin SVD | Con SVD |
|---------|---------|---------|
| **Número de redes** | 8 | k (10-20) |
| **Input por red** | `[r, x, y]` (3D) | `[r]` (1D) |
| **Output por red** | `[Tx_i, Ty_i]` (2D por nodo) | `V(r,i)` (escalar) |
| **Parámetros de salida** | `83,232` | `k ≈ 10-20` |
| **Complejidad** | Alta | Baja |
| **Precisión** | Exacta | Aproximada (controlable) |
| **Velocidad de predicción** | Lenta | Rápida |

## 8. Código Detallado - Línea por Línea

### Líneas 6-20: Construcción de T_SVD

```matlab
r = 0:0.02:0.999;  % 50 radios espaciados uniformemente

for i=1:size(r,2)  % Loop sobre 50 radios
    % Construir nombre de archivo
    string = strrep("UL_r"+num2str(r(i), '%.4f'), ".", "_")+"-50x50"+".mat";
    FileName = fullfile('NilTFMFiles','DataVariables','50x50',string);
    
    % Cargar T y mesh para este radio
    load(FileName,"T","mesh");
    
    % Pre-allocar T_SVD en la primera iteración
    if i==1
        % T_SVD: [nDofs×8 filas × nRadios columnas]
        T_SVD = zeros(mesh.nnodes*mesh.ndim*8, size(r,2));
    end
    
    % Vectorizar T y almacenar en columna i
    T_SVD(:,i) = T(:);  % T(:) convierte [nDofs×8] → [nDofs×8×1]
end
```

**Resultado**: `T_SVD` es `[41,616 × 50]`

### Línea 22: Descomposición SVD

```matlab
[U,S,V] = svd(T_SVD,'econ');
```

**`'econ'`**: Modo económico (solo calcula las primeras `min(m,n)` columnas de `U` y `V`)

**Resultado**:
- `U`: `[41,616 × 50]` (modos espaciales)
- `S`: `[50 × 50]` (valores singulares, diagonal)
- `V`: `[50 × 50]` (modos paramétricos)

### Líneas 27-44: Visualización de V (agrupada)

```matlab
step = 10;  % Mostrar 10 modos por ventana
Nwindow = ceil(size(V,2)/step);  % Número de ventanas

for j=1:Nwindow  % Para cada ventana
    figure('Position',[75 100 1400 600]);
    tiledlayout(2,5);  % Grid 2×5
    
    for i=1:step  % 10 subplots por ventana
        ax = nexttile;
        plot(r, V(:,idx), 'LineWidth', 1.5);  % V(:,idx) vs r
        xlabel('r');
        ylabel("V(:,"+idx);
        title("V-"+ idx);
        grid on
        idx = idx+1;
    end
end
```

**Muestra**: Cómo varían los primeros 50 modos paramétricos con el radio

### Líneas 46-55: Visualización de V (todos)

```matlab
tiledlayout(10,5);  % Grid 10×5 = 50 subplots
for i=1:50
    ax = nexttile;
    plot(r, V(:,i), 'LineWidth', 1.5);
    xlabel('r');
    ylabel("V"+i);
    title("V"+ i);
    grid on
    i = i+1;  % Nota: esto es redundante (i se incrementa en el for)
end
```

**Muestra**: Todos los 50 modos paramétricos en una sola figura

### Líneas 57-62: Visualización de S

```matlab
figure
plot(log(diag(S)), 'LineWidth', 1.5);
title("S singular values");
ylabel('value');
xlabel('column');
```

**Muestra**: Decaimiento de los valores singulares (en escala logarítmica)

**Interpretación**:
- Si decae rápido → pocos modos dominantes
- Si decae lento → muchos modos necesarios

### Líneas 64-74: Visualización de U (comentado)

```matlab
% Código comentado que visualizaría U como campos sobre la malla
% s.mesh = mesh;
% s.order = 'P1';
% n = 10;
% 
% Ufun = cell(1,n);
% for i=1:n
%     s.fValues = U(:,i);  % Modo espacial i
%     Ufun{1,i} = LagrangianFunction(s);
%     Ufun{1,i}.plot();  % Visualiza como campo sobre la malla
% end
```

**Si se descomentara**: Mostraría los primeros 10 modos espaciales como campos de desplazamiento sobre la malla

### Líneas 76-83: Visualización de U (valores)

```matlab
figure('Position',[75 100 1400 600]);
tiledlayout(2,5);  % Grid 2×5
for i=1:10
    ax = nexttile;
    plot(U(:,i));  % Valores del modo espacial i
    ylabel("U values");
    title("U column "+ i);
end
```

**Muestra**: Los primeros 10 modos espaciales como vectores (índice vs valor)

## 9. Aplicaciones Prácticas

### 9.1. Selección del Número de Modos

**Criterio basado en energía**:
```matlab
energy_ratio = cumsum(diag(S).^2) / sum(diag(S).^2);
k = find(energy_ratio > 0.99, 1);  % 99% de energía
```

**Criterio basado en tolerancia**:
```matlab
tol = 1e-6;
k = find(diag(S) / S(1,1) > tol, 1, 'last');
```

### 9.2. Reconstrucción para un Radio Nuevo

```matlab
% Predicción de V para r_new usando red neuronal
V_predicted = predict_V_from_NN(r_new);  % [k × 1]

% Reconstrucción de T
T_reconstructed = U(:,1:k) * S(1:k,1:k) * V_predicted;
T_reconstructed = reshape(T_reconstructed, [nDofs, 8]);
```

### 9.3. Compresión de Datos

**Almacenamiento original**:
- 50 archivos `.mat` con `T` completo
- Tamaño: `50 × (nDofs × 8 × 8 bytes) ≈ 16.6 MB`

**Almacenamiento con SVD**:
- `U(:,1:k)`: `[nDofs×8 × k]`
- `S(1:k,1:k)`: `[k × k]`
- `V`: `[50 × k]`
- Tamaño: `(nDofs×8 × k + k² + 50×k) × 8 bytes`
- Para `k=20`: `≈ 6.7 MB` (reducción de ~60%)

## 10. Ventajas y Limitaciones

### Ventajas

1. **Reducción masiva de dimensionalidad**: De 83,232 → 10-20 parámetros
2. **Velocidad**: Predicción mucho más rápida
3. **Robustez**: Modos dominantes capturan la física esencial
4. **Interpretabilidad**: Modos tienen significado físico claro

### Limitaciones

1. **Aproximación**: No es exacta (pero controlable)
2. **Dependencia de datos**: SVD depende de los radios de entrenamiento
3. **Extrapolación**: Puede fallar fuera del rango de `r` usado
4. **Linealidad**: Asume que la variación es principalmente lineal en el espacio de modos

## 11. Relación con Option 3

**Option 3 sin SVD**:
- Entrenar 8 redes: `[r, x, y] → [Tx_i, Ty_i]`
- Predicción directa de desplazamientos

**Option 3 con SVD**:
- Entrenar k redes: `r → V(r,i)`
- Reconstrucción: `T = U × S × V'`
- Predicción indirecta pero más eficiente

**Elección**:
- **Sin SVD**: Mayor precisión, más lento, más datos necesarios
- **Con SVD**: Menor precisión (controlable), más rápido, menos datos necesarios


