# Análisis en Profundidad: Generación de Datos para Entrenamiento de Redes Neuronales

## 1. Visión General del Proceso

El proceso genera datos de entrenamiento para dos tipos de redes neuronales:
- **Redes para T (Option 3)**: Predicen la matriz de downscaling `T` que mapea desplazamientos del espacio grueso (8 modos) al espacio fino
- **Redes para K (Option 2)**: Predicen la matriz de rigidez gruesa `Kcoarse` (8×8)

## 2. Flujo Detallado del Proceso Actual

### 2.1. DataGeneration_Abril.m - Script Principal

```matlab
r = 0:0.1:0.999;  % Vector de radios: [0, 0.1, 0.2, ..., 0.9]
nelem = 40;       % Elementos por lado de la malla
```

**Para cada radio `r(j)`:**
1. Llama a `LevelSetInclusionAuto_abril(r(j), 1, nelem, doplot)`
2. Recibe:
   - `u`: Matriz T (desplazamientos para 8 modos) - tamaño: `[nDofs × 8]`
   - `L`: Multiplicadores de Lagrange
   - `mesh`: Malla del dominio
   - `Kcoarse`: Matriz de rigidez gruesa `[8 × 8]`

3. **Procesamiento de T (u)**:
   - Cada columna de `u` representa un modo: `u(:,j)` = modo j
   - Cada modo tiene componentes x e y: `u = [u_x1, u_y1, u_x2, u_y2, ..., u_x8, u_y8]`
   - Se reshape para agrupar (Tx, Ty) por nodo:
     ```matlab
     t1 = reshape(u(:,1).', 2, []).';  % [nnodes × 2] con [Tx1, Ty1]
     ```
   - Se crea formato CSV: `[r, x, y, Tx1, Ty1, Tx2, Ty2, ..., Tx8, Ty8]`
   - Cada fila = un nodo con su radio, coordenadas y desplazamientos de los 8 modos

4. **Procesamiento de K**:
   - `Kcoarse` es simétrica 8×8 → 36 componentes únicas
   - Se extrae triangular superior y se vectoriza
   - Formato CSV: `[r, K11, K12, K13, ..., K88]` (36 componentes)

### 2.2. LevelSetInclusionAuto_abril.m - Solución del Problema Elástico

#### 2.2.1. Inicialización
```matlab
obj.init(r, i, nelem)
```
- `r`: Radio de la inclusión
- `i`: Dirección del nodo (1-8) para condiciones de contorno
- `nelem`: Número de elementos por lado

#### 2.2.2. Creación de Malla
```matlab
createMesh() → createReferenceMesh()
```
- Crea malla estructurada en `[-1, 1] × [-1, 1]`
- **NOTA**: Actualmente NO crea agujero físico (líneas comentadas)
- Malla: `nelem × nelem` elementos triangulares

#### 2.2.3. Material con Inclusión
```matlab
computeElasticProperties()
```
- **E1 = 1** (material exterior, rígido)
- **E2 = E1/1000 = 0.001** (material interior, blando)
- **nu = 1/3** (Poisson constante)
- Función analítica:
  ```matlab
  E(x,y) = (dist < r) * E2 + (dist >= r) * E1
  ```
- La inclusión es un círculo centrado en (0,0) con radio `r`

#### 2.2.4. Condiciones de Contorno para 8 Modos Gruesos

**Física**: Se resuelven 8 problemas elásticos independientes, uno por cada modo grueso.

**Modos Gruesos (8 modos)**:
- **Modo 1**: Desplazamiento horizontal en esquina 1
- **Modo 2**: Desplazamiento vertical en esquina 1
- **Modo 3**: Desplazamiento horizontal en esquina 2
- **Modo 4**: Desplazamiento vertical en esquina 2
- **Modo 5**: Desplazamiento horizontal en esquina 3
- **Modo 6**: Desplazamiento vertical en esquina 3
- **Modo 7**: Desplazamiento horizontal en esquina 4
- **Modo 8**: Desplazamiento vertical en esquina 4

**Implementación**:
```matlab
createBoundaryConditions()
```
- `nodeDirection = i` (1-8) determina qué modo se resuelve
- **4 esquinas** del dominio: (xMin,yMin), (xMax,yMin), (xMin,yMax), (xMax,yMax)
- **Para cada modo**:
  - 3 esquinas: desplazamiento = 0 (Dirichlet)
  - 1 esquina: desplazamiento = 1 en dirección específica (Dirichlet)
  - Sin fuerzas externas (RHS = 0)

**Matriz de asignación**:
```matlab
assignMatrix = [2 1 0 0 0 0 0 0  % Modo 1: esquina 1, dirección x
                0 0 2 1 0 0 0 0  % Modo 2: esquina 1, dirección y
                0 0 0 0 2 1 0 0  % Modo 3: esquina 2, dirección x
                0 0 0 0 0 0 2 1  % Modo 4: esquina 2, dirección y
                1 2 1 2 1 2 1 2];% Última fila: dirección del desplazamiento unitario
```

#### 2.2.5. Solución del Problema Elástico

**Sistema de ecuaciones**:
```
[K  C] [u]   [0]
[C' 0] [L] = [rdir]
```

Donde:
- **K**: Matriz de rigidez fina `[nDofs × nDofs]`
- **C**: Matriz de restricciones (Lagrange) `[nDofs × nConstraints]`
- **u**: Desplazamientos `[nDofs × 8]` (8 columnas = 8 modos)
- **L**: Multiplicadores de Lagrange
- **rdir**: RHS para condiciones Dirichlet (8 funciones base)

**Cálculo de K**:
```matlab
K = IntegrateLHS(@(u,v) DDP(SymGrad(v), DDP(C, SymGrad(u))), ...)
```
- Integración de forma débil: `∫ ε(v) : C : ε(u) dΩ`
- `C`: Tensor de material (Young's modulus variable)

**Cálculo de C (restricciones)**:
```matlab
C = IntegrateLHS(@(u,v) DP(v, u), test, dLambda, mesh, 'Boundary', 2)
```
- Restricciones en el contorno: `∫ v · u dΓ`

**RHS Dirichlet (rdir)**:
```matlab
RHSdirichlet()
```
- 8 funciones base bilineales en las esquinas:
  - `f1x, f1y`: Esquina 1 (xMin, yMin)
  - `f2x, f2y`: Esquina 2 (xMax, yMin)
  - `f3x, f3y`: Esquina 3 (xMax, yMax)
  - `f4x, f4y`: Esquina 4 (xMin, yMax)
- Cada función es 1/4 del producto de funciones de forma bilineales

**Solución**:
```matlab
sol = LHS \ RHS
u = sol(1:nDofs, :)  % 8 columnas, una por modo
```

#### 2.2.6. Cálculo de Kcoarse

```matlab
Kcoarse = u.' * obj.stiffness * u
```

**Interpretación física**:
- `Kcoarse[i,j] = u_i^T * K_fine * u_j`
- Representa la rigidez del sistema proyectada al espacio grueso
- `u_i`: i-ésima columna de T (modo i)
- `Kcoarse`: Matriz 8×8 simétrica

## 3. Estructura de Datos Generados

### 3.1. DataT.csv (para entrenar redes de T)

**Formato**: `[r, x, y, Tx1, Ty1, Tx2, Ty2, ..., Tx8, Ty8]`

**Dimensiones**:
- Filas: `nRadii × nnodes` (ej: 10 × 1681 = 16,810 filas para 40×40 elementos)
- Columnas: 19 (1 radio + 2 coordenadas + 16 componentes de desplazamiento)

**Uso**: Entrenar 8 redes neuronales separadas:
- Red 1: Input `[r, x, y]` → Output `[Tx1, Ty1]`
- Red 2: Input `[r, x, y]` → Output `[Tx2, Ty2]`
- ...
- Red 8: Input `[r, x, y]` → Output `[Tx8, Ty8]`

### 3.2. DataK.csv (para entrenar red de K)

**Formato**: `[r, K11, K12, K13, ..., K88]` (36 componentes únicas)

**Dimensiones**:
- Filas: `nRadii` (ej: 10 filas)
- Columnas: 37 (1 radio + 36 componentes de K)

**Uso**: Entrenar 1 red neuronal:
- Input: `[r]`
- Output: `[K11, K12, ..., K88]` (36 valores)

## 4. Problemas y Limitaciones del Código Actual

1. **Malla sin agujero físico**: La inclusión es solo una propiedad de material, no un agujero geométrico
2. **Proceso secuencial**: Resuelve un radio a la vez (no paralelizable fácilmente)
3. **Almacenamiento redundante**: Guarda .mat individuales + CSV
4. **Código no modular**: Todo mezclado en una clase
5. **Sin validación**: No verifica que los datos sean consistentes

## 5. Propuesta de Estructura Optimizada

### 5.1. Arquitectura Modular Propuesta

```
DataGeneration/
├── Core/
│   ├── RVESolver.m              % Soluciona problema elástico para un RVE
│   ├── MaterialDefinition.m     % Define propiedades de material con inclusión
│   ├── BoundaryConditions.m     % Crea condiciones de contorno para modos gruesos
│   └── CoarseModeSolver.m       % Resuelve los 8 modos gruesos
├── DataProcessing/
│   ├── TDataProcessor.m         % Procesa y formatea datos de T
│   ├── KDataProcessor.m         % Procesa y formatea datos de K
│   └── DataValidator.m          % Valida consistencia de datos
├── Export/
│   ├── CSVExporter.m            % Exporta a CSV
│   └── MATExporter.m            % Exporta a .mat (opcional)
└── DataGenerationPipeline.m    % Script principal orquestador
```

### 5.2. Flujo Optimizado Propuesto

1. **Preparación**:
   - Definir rango de radios
   - Crear malla de referencia (una vez)
   - Pre-allocar arrays

2. **Loop paralelizable** (parfor):
   - Para cada radio:
     - Crear material con inclusión
     - Resolver 8 modos gruesos (puede ser paralelo)
     - Calcular Kcoarse
     - Almacenar resultados

3. **Post-procesamiento**:
   - Formatear datos de T
   - Formatear datos de K
   - Validar datos
   - Exportar a CSV

### 5.3. Optimizaciones Propuestas

1. **Paralelización**: Usar `parfor` para diferentes radios
2. **Reutilización de malla**: Crear malla una vez, reutilizar
3. **Almacenamiento eficiente**: Solo CSV, eliminar .mat individuales
4. **Validación automática**: Verificar simetría de K, ortogonalidad de modos
5. **Caché de resultados**: Guardar resultados intermedios para evitar recálculo

## 6. Implementación: OptimizedDataGeneration.m

### 6.1. Arquitectura de la Clase

La clase `OptimizedDataGeneration` implementa todas las optimizaciones propuestas:

**Estructura**:
- **Constructor**: Inicializa parámetros y crea malla de referencia (una vez)
- **generateData()**: Orquesta el proceso completo con paralelización opcional
- **solveForRadius()**: Resuelve problema elástico para un radio (llamado en loop)
- **solveCoarseModes()**: Resuelve los 8 modos gruesos
- **processTData()**: Formatea datos de T para CSV
- **processKData()**: Formatea datos de K para CSV
- **exportToCSV()**: Exporta a archivos CSV
- **validateData()**: Valida consistencia de datos

### 6.2. Comparación con Código Original

| Aspecto | DataGeneration_Abril.m | OptimizedDataGeneration.m |
|---------|------------------------|---------------------------|
| **Paralelización** | No | Sí (parfor opcional) |
| **Reutilización de malla** | No (crea cada vez) | Sí (una vez en constructor) |
| **Almacenamiento** | .mat individuales + CSV | Solo CSV (opcional .mat) |
| **Validación** | Manual | Automática |
| **Modularidad** | Baja (todo en una clase) | Alta (métodos separados) |
| **Agujeros físicos** | No (comentado) | Opcional (parámetro) |
| **Código** | ~94 líneas script | ~450 líneas clase modular |

### 6.3. Uso de la Clase Optimizada

```matlab
% Configuración
radii = 0:0.1:0.9;  % 10 radios
nelem = 40;          % 40×40 elementos

% Crear generador (sin agujeros físicos - más rápido)
dataGen = OptimizedDataGeneration(radii, nelem, 'usePhysicalHoles', false);

% Generar datos (paralelo si está disponible)
dataGen.generateData('parallel', true);

% Validar datos
dataGen.validateData();

% Exportar a CSV
dataGen.exportToCSV('AbrilTFGfiles');
```

### 6.4. Detalles de Implementación Clave

#### 6.4.1. Resolución de Modos Gruesos

**Problema**: Resolver 8 problemas elásticos independientes:
```
[K_fine  C] [u_i]   [0     ]
[C'      0] [L_i] = [rdir_i]
```

**Estrategia**:
1. Calcular `K_fine` una vez (reutilizable para todos los modos)
2. Para cada modo i (1-8):
   - Crear BC específica: 3 esquinas con u=0, 1 esquina con u=1
   - Calcular `C`: matriz de restricciones (Lagrange) en el contorno
   - Calcular `rdir_i`: RHS Dirichlet usando función base bilineal i
   - Resolver sistema aumentado
   - Extraer `u_i` (columna i de la matriz T)

#### 6.4.2. Funciones Base Bilineales

Las 8 funciones base representan desplazamientos unitarios en las esquinas:

- **Modo 1**: Esquina 1 (xMin, yMin), dirección x → `f1x`
- **Modo 2**: Esquina 1 (xMin, yMin), dirección y → `f1y`
- **Modo 3**: Esquina 2 (xMax, yMin), dirección x → `f2x`
- **Modo 4**: Esquina 2 (xMax, yMin), dirección y → `f2y`
- **Modo 5**: Esquina 3 (xMax, yMax), dirección x → `f3x`
- **Modo 6**: Esquina 3 (xMax, yMax), dirección y → `f3y`
- **Modo 7**: Esquina 4 (xMin, yMax), dirección x → `f4x`
- **Modo 8**: Esquina 4 (xMin, yMax), dirección y → `f4y`

**Forma matemática**:
```
f1x(x,y) = [1/4 * (1-x) * (1-y), 0]
f1y(x,y) = [0, 1/4 * (1-x) * (1-y)]
f2x(x,y) = [1/4 * (1+x) * (1-y), 0]
...
```

Estas funciones son productos de funciones de forma bilineales estándar.

#### 6.4.3. Formato de Datos T

**Estructura original de u**:
- `u`: `[nDofs × 8]` donde `nDofs = nnodes × ndim = nnodes × 2`
- Cada columna i: `[u_x1, u_y1, u_x2, u_y2, ..., u_x_nnodes, u_y_nnodes]'`

**Transformación a formato CSV**:
1. Para cada modo i:
   - Reshape: `reshape(u(:,i), 2, nnodes)'` → `[nnodes × 2]`
   - Resultado: `[Tx_i, Ty_i]` por nodo
2. Concatenar todos los modos: `[Tx1, Ty1, Tx2, Ty2, ..., Tx8, Ty8]`
3. Agregar radio y coordenadas: `[r, x, y, Tx1, Ty1, ..., Tx8, Ty8]`
4. Cada fila = un nodo con sus desplazamientos para los 8 modos

**Uso para entrenamiento**:
- Entrenar 8 redes neuronales separadas
- Red i: Input `[r, x, y]` → Output `[Tx_i, Ty_i]`

#### 6.4.4. Formato de Datos K

**Kcoarse**: Matriz simétrica 8×8 → 36 componentes únicas

**Extracción**:
1. Obtener triangular superior: `triu(K)`
2. Vectorizar: `[K11, K12, K13, ..., K18, K22, K23, ..., K88]`
3. Formato CSV: `[r, K11, K12, ..., K88]` (37 columnas)

**Uso para entrenamiento**:
- Entrenar 1 red neuronal
- Input: `[r]`
- Output: `[K11, K12, ..., K88]` (36 valores)

### 6.5. Ventajas de la Implementación Optimizada

1. **Rendimiento**:
   - Paralelización reduce tiempo total
   - Reutilización de malla evita recálculos
   - Pre-allocación de arrays mejora eficiencia de memoria

2. **Mantenibilidad**:
   - Código modular y bien documentado
   - Métodos pequeños y enfocados
   - Fácil de extender y modificar

3. **Robustez**:
   - Validación automática de datos
   - Manejo de errores mejorado
   - Opciones configurables

4. **Flexibilidad**:
   - Soporte opcional para agujeros físicos
   - Fácil cambiar parámetros
   - Exportación configurable

