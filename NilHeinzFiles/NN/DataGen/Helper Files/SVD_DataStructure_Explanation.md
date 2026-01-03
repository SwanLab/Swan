# Estructura de Datos después de SVD: Explicación Detallada

## 1. ¿Qué forma tiene T después de aplicar SVD?

**Respuesta corta**: T **NO se transforma directamente**. En su lugar, se **descompone** en tres componentes (U, S, V) y se **reconstruye** cuando se necesita.

### Antes del SVD:
- **T original**: Matriz `[nDofs × 8]` para cada radio `r`
  - `nDofs = nnodes × 2` (grados de libertad)
  - 8 columnas = 8 modos gruesos
  - Ejemplo: Para malla 50×50, `T` es `[5202 × 8]`

### Después del SVD:
- **T se descompone en**: `T_SVD = U × S × V'`
- **T se reconstruye cuando se necesita**: `T_reconstructed = U × S × V_predicted'`

## 2. Estructura de U, S, V

### U: Modos Espaciales (NO es una función, es una matriz)

**Forma**: `U` es una matriz `[nDofs×8 × k]` donde:
- `nDofs×8` = número total de componentes de T (41,616 para malla 50×50)
- `k` = número de modos retenidos (típicamente 10-20)

**Interpretación**:
- Cada **columna** `U(:,i)` = un modo espacial (patrón de desplazamiento fijo en el espacio)
- `U` es **fijo** (pre-computado una vez durante el SVD)
- `U` NO depende del radio `r`

**Ejemplo**:
```matlab
% U tiene forma [41616 × k]
% U(:,1) = primer modo espacial (vector de 41616 componentes)
% U(:,2) = segundo modo espacial
% ...
% U(:,k) = k-ésimo modo espacial
```

**Inputs/Outputs de U**:
- **Input**: Índice del modo `i` (1 a k)
- **Output**: Vector `U(:,i)` de tamaño `[nDofs×8 × 1]`
- **U NO es una función**, es una **matriz pre-computada**

### S: Valores Singulares (matriz diagonal)

**Forma**: `S` es una matriz diagonal `[k × k]`

**Interpretación**:
- `S(i,i) = σᵢ` = valor singular del modo `i`
- `σ₁ ≥ σ₂ ≥ ... ≥ σₖ`
- `S` es **fijo** (pre-computado una vez durante el SVD)
- `S` NO depende del radio `r`

**Ejemplo**:
```matlab
% S es [k × k] diagonal
% S = diag([σ₁, σ₂, ..., σₖ])
```

### V: Modos Paramétricos (se predice con NN)

**Forma**: `V` es una matriz `[nRadii × k]` donde:
- `nRadii` = número de radios en el conjunto de entrenamiento (ej: 50)
- `k` = número de modos retenidos

**Interpretación**:
- Cada **fila** `V(j,:)` = coeficientes de los k modos para el radio `r(j)`
- Cada **columna** `V(:,i)` = cómo varía el modo `i` con el radio `r`
- `V` se **predice** con redes neuronales para radios nuevos

**Ejemplo**:
```matlab
% V tiene forma [50 × k] (para 50 radios de entrenamiento)
% V(1,:) = [V₁(r₁), V₂(r₁), ..., Vₖ(r₁)]  % Para radio r₁
% V(2,:) = [V₁(r₂), V₂(r₂), ..., Vₖ(r₂)]  % Para radio r₂
% ...
% V(:,1) = [V₁(r₁), V₁(r₂), ..., V₁(r₅₀)]'  % Modo 1 para todos los radios
```

## 3. Flujo de Predicción con NN

### Paso 1: Entrenamiento de Redes Neuronales

**Entrenar k redes neuronales** (una por cada modo):

```matlab
% Para cada modo i = 1, 2, ..., k:
%   Red NN_i: Input = r (escalar)
%            Output = V(r, i) (escalar)
```

**Datos de entrenamiento**:
- Input: `r` (vector de radios: `[r₁, r₂, ..., r₅₀]`)
- Output: `V(:,i)` (columna i de V: `[V(r₁,i), V(r₂,i), ..., V(r₅₀,i)]'`)

**Formato CSV para entrenamiento**:
```
r, V1, V2, ..., Vk
r₁, V(r₁,1), V(r₁,2), ..., V(r₁,k)
r₂, V(r₂,1), V(r₂,2), ..., V(r₂,k)
...
r₅₀, V(r₅₀,1), V(r₅₀,2), ..., V(r₅₀,k)
```

### Paso 2: Predicción para un Radio Nuevo

```matlab
% Dado un radio nuevo r_new:

% 1. Predecir V usando las k redes neuronales
V_predicted = zeros(k, 1);
for i = 1:k
    V_predicted(i) = NN_i.predict(r_new);  % Predicción del modo i
end
% V_predicted es [k × 1]

% 2. Reconstruir T
T_vec = U(:,1:k) * S(1:k,1:k) * V_predicted;  % [nDofs×8 × 1]

% 3. Reshape a forma original
T_reconstructed = reshape(T_vec, [nDofs, 8]);  % [nDofs × 8]
```

## 4. Comparación: U vs Función

### U NO es una función

**U es una matriz pre-computada**:
- Se calcula **una vez** durante el SVD
- Se almacena como matriz `[nDofs×8 × k]`
- **NO necesita inputs** (ya está calculado)
- Se **usa directamente** en la reconstrucción

**Ejemplo de uso**:
```matlab
% U ya está calculado y almacenado
% No necesitas llamar una función, solo indexar:
U_mode_1 = U(:, 1);  % Primer modo espacial
U_mode_2 = U(:, 2);   % Segundo modo espacial
```

### V SÍ se predice con funciones (NN)

**V se predice con redes neuronales**:
- Se entrena **k redes neuronales** (una por modo)
- Cada red es una **función**: `V_i(r) = NN_i(r)`
- Input: `r` (escalar)
- Output: `V(r,i)` (escalar)

**Ejemplo de uso**:
```matlab
% Para un radio nuevo r_new:
V_predicted = zeros(k, 1);
for i = 1:k
    V_predicted(i) = NN_i.computeOutputValues(r_new);
end
```

## 5. Estructura Completa de Datos


**Después de ejecutar, tienes**:

1. **U**: Matriz `[nDofs×8 × k]` - Modos espaciales (fijos)
2. **S**: Matriz diagonal `[k × k]` - Valores singulares (fijos)
3. **V**: Matriz `[nRadii × k]` - Modos paramétricos (para entrenamiento)

### Archivos para Entrenamiento de NN

**VTrainingData.csv** (generado por TrainingDataGenerator):
```
r, V1, V2, ..., Vk
0.00, V(0.00,1), V(0.00,2), ..., V(0.00,k)
0.02, V(0.02,1), V(0.02,2), ..., V(0.02,k)
...
0.98, V(0.98,1), V(0.98,2), ..., V(0.98,k)
```

**Formato**:
- Columna 1: `r` (input para NN)
- Columnas 2 a k+1: `V1, V2, ..., Vk` (outputs para entrenar k redes)

### Archivos para Reconstrucción

**SVD_Results.mat** (generado por TrainingDataGenerator):
```matlab
% Contiene:
U        % [nDofs×8 × k] - Modos espaciales
S        % [k × k] - Valores singulares (diagonal)
V        % [nRadii × k] - Modos paramétricos (solo para referencia)
k        % Número de modos retenidos
radii    % Radios de entrenamiento
meshRef  % Malla de referencia (opcional)
```

## 6. Ejemplo Práctico de Reconstrucción

```matlab
% 1. Cargar resultados SVD
load('SVD_Results.mat', 'U', 'S', 'k', 'meshRef');

% 2. Cargar redes neuronales entrenadas
load('NN_V1.mat', 'opt1');  % Red para modo 1
load('NN_V2.mat', 'opt2');  % Red para modo 2
% ... (k redes en total)

% 3. Dado un radio nuevo r_new = 0.35
r_new = 0.35;

% 4. Predecir V usando las k redes
V_predicted = zeros(k, 1);
V_predicted(1) = opt1.computeOutputValues(r_new);
V_predicted(2) = opt2.computeOutputValues(r_new);
% ... (para los k modos)

% 5. Reconstruir T
T_vec = U * S * V_predicted;  % [nDofs×8 × 1]

% 6. Reshape a forma original
nDofs = size(meshRef.coord, 1) * meshRef.ndim;
T_reconstructed = reshape(T_vec, [nDofs, 8]);  % [nDofs × 8]

% 7. T_reconstructed es ahora la matriz T para el radio r_new
%    Cada columna T_reconstructed(:,i) es el desplazamiento fino para el modo grueso i
```

## 7. Resumen: Inputs y Outputs

### U (Modos Espaciales)
- **Tipo**: Matriz pre-computada
- **Forma**: `[nDofs×8 × k]`
- **Input**: Índice del modo `i` (1 a k)
- **Output**: Vector `U(:,i)` de tamaño `[nDofs×8 × 1]`
- **Dependencia**: **NO depende de r** (fijo)

### S (Valores Singulares)
- **Tipo**: Matriz diagonal pre-computada
- **Forma**: `[k × k]`
- **Input**: Ninguno (se usa directamente)
- **Output**: Matriz diagonal con valores singulares
- **Dependencia**: **NO depende de r** (fijo)

### V (Modos Paramétricos)
- **Tipo**: Se predice con redes neuronales
- **Forma**: `[k × 1]` para un radio dado
- **Input**: Radio `r` (escalar)
- **Output**: Vector `V(r) = [V(r,1), V(r,2), ..., V(r,k)]'`
- **Dependencia**: **SÍ depende de r** (se predice)

### T (Reconstruido)
- **Tipo**: Matriz reconstruida
- **Forma**: `[nDofs × 8]`
- **Input**: Radio `r` (para predecir V)
- **Output**: Matriz T completa
- **Reconstrucción**: `T = reshape(U × S × V(r), [nDofs, 8])`

## 8. Ventajas de esta Estructura

1. **U y S son fijos**: Se calculan una vez y se reutilizan
2. **Solo V se predice**: Reducción masiva de dimensionalidad (de 83,232 → k ≈ 10-20)
3. **Reconstrucción rápida**: Multiplicación de matrices es muy eficiente
4. **Flexibilidad**: Puedes cambiar k (número de modos) sin regenerar U y S

