# Red Neuronal para Predicción de Matrices de Rigidez 8x8

Este tutorial muestra cómo entrenar y usar una red neuronal para predecir matrices de rigidez 8x8 basadas en un parámetro de entrada único.

## Estructura del Sistema

### Arquitectura de la Red Neuronal
- **Entrada**: 1 característica (parámetro de entrada, ej: densidad, espesor, etc.)
- **Salida**: 64 características (matriz de rigidez 8x8 aplanada)
- **Arquitectura**: Red profunda con 8 capas ocultas de 256 neuronas cada una

### Archivos Principales

1. **`trainNN_stiffness_matrix.m`** - Script principal de entrenamiento
2. **`StiffnessMatrixHandler.m`** - Utilidades para manejo de matrices de rigidez
3. **`genDB_stiffness_matrix.m`** - Generador de datos de entrenamiento
4. **`OptimizeStiffnessMatrixTutorial.m`** - Optimización usando la red entrenada

## Uso del Sistema

### 1. Generar Datos de Entrenamiento

```matlab
% Ejecutar el generador de datos
genDB_stiffness_matrix();
```

Esto creará un archivo CSV con:
- Columna 1: parámetro de entrada
- Columnas 2-65: elementos de la matriz 8x8 aplanada

### 2. Entrenar la Red Neuronal

```matlab
% Ejecutar el script de entrenamiento
trainNN_stiffness_matrix();
```

### 3. Usar la Red para Predicciones

```matlab
% Cargar la red entrenada
load('Tutorials/ChomogNetworkTutorial/Networks/network_stiffness_matrix.mat');

% Hacer una predicción
input_param = 1.5;
K_vector = opt.computeOutputValues(input_param);
K_matrix = StiffnessMatrixHandler.vectorToMatrix(K_vector);

% Visualizar la matriz
StiffnessMatrixHandler.plotStiffnessMatrix(K_matrix);
```

### 4. Optimización

```matlab
% Ejecutar optimización
OptimizeStiffnessMatrixTutorial();
```

## Funciones de Utilidad

### StiffnessMatrixHandler

#### Conversión de Formatos
```matlab
% Vector a matriz
K = StiffnessMatrixHandler.vectorToMatrix(vector_64);

% Matriz a vector
vector = StiffnessMatrixHandler.matrixToVector(K_8x8);
```

#### Validación de Matrices
```matlab
% Verificar si es una matriz de rigidez válida
isValid = StiffnessMatrixHandler.validateStiffnessMatrix(K);

% Hacer que una matriz sea válida (simétrica y definida positiva)
K_valid = StiffnessMatrixHandler.makeValidStiffnessMatrix(K);
```

#### Análisis
```matlab
% Analizar propiedades de la matriz
[condition_number, max_eigenval, min_eigenval] = StiffnessMatrixHandler.analyzeMatrix(K);

% Visualizar matriz
StiffnessMatrixHandler.plotStiffnessMatrix(K, 'Título');
```

## Tipos de Optimización Disponibles

1. **`maxStiffness`** - Maximizar la traza de la matriz
2. **`minStiffness`** - Minimizar la traza de la matriz
3. **`maxConditionNumber`** - Maximizar el número de condición
4. **`minConditionNumber`** - Minimizar el número de condición
5. **`targetStiffness`** - Acercarse a una matriz objetivo
6. **`maxDiagonal`** - Maximizar elementos diagonales
7. **`minDiagonal`** - Minimizar elementos diagonales

## Personalización

### Modificar la Función de Generación de Datos

En `genDB_stiffness_matrix.m`, modifica la función `generateStiffnessMatrix()` para usar tu modelo físico específico:

```matlab
function K = generateStiffnessMatrix(input_param)
    % Tu modelo físico aquí
    % Ejemplo: K = tu_modelo_fisico(input_param);
end
```

### Ajustar Hiperparámetros

En `trainNN_stiffness_matrix.m`:

```matlab
% Ajustar arquitectura de la red
hiddenLayers = 256 .* ones(1, 8);  % 8 capas de 256 neuronas

% Ajustar parámetros de entrenamiento
learningRate = 0.001;
maxEpochs = 1000;
```

### Agregar Nuevos Tipos de Optimización

En `OptimizeStiffnessMatrixTutorial.m`, agrega nuevos casos en la función `stiffnessMatrixCost()`:

```matlab
case 'nuevoTipo'
    % Tu función de costo aquí
    J = tu_funcion_costo(K);
    grad = tu_gradiente(dK);
```

## Consideraciones Importantes

1. **Simetría**: Las matrices de rigidez deben ser simétricas
2. **Definida Positiva**: Deben tener todos los eigenvalores positivos
3. **Escalabilidad**: Para matrices más grandes, considera usar representaciones más eficientes
4. **Validación**: Siempre valida que las matrices generadas sean físicamente válidas

## Ejemplo Completo

```matlab
% 1. Generar datos
genDB_stiffness_matrix();

% 2. Entrenar red
trainNN_stiffness_matrix();

% 3. Cargar red entrenada
load('Tutorials/ChomogNetworkTutorial/Networks/network_stiffness_matrix.mat');

% 4. Hacer predicción
input_param = 1.2;
K_vector = opt.computeOutputValues(input_param);
K_matrix = StiffnessMatrixHandler.vectorToMatrix(K_vector);

% 5. Validar y visualizar
isValid = StiffnessMatrixHandler.validateStiffnessMatrix(K_matrix);
fprintf('Matriz válida: %s\n', mat2str(isValid));

StiffnessMatrixHandler.plotStiffnessMatrix(K_matrix, 'Matriz Predicha');

% 6. Optimizar
OptimizeStiffnessMatrixTutorial();
```

## Troubleshooting

### Problemas Comunes

1. **Matrices no simétricas**: Usa `StiffnessMatrixHandler.enforceSymmetry()`
2. **Matrices no definidas positivas**: Usa `StiffnessMatrixHandler.enforcePositiveDefinite()`
3. **Convergencia lenta**: Ajusta `learningRate` y `maxEpochs`
4. **Overfitting**: Aumenta `lambda` para regularización L2

### Validación de Resultados

```matlab
% Verificar propiedades físicas
eigenvals = eig(K_matrix);
fprintf('Eigenvalores positivos: %s\n', mat2str(all(eigenvals > 0)));
fprintf('Simétrica: %s\n', mat2str(issymmetric(K_matrix)));
fprintf('Número de condición: %.2e\n', cond(K_matrix));
```
