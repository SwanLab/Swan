# Celdas biestables: Chen 2021 y Rafsanjani & Pasini 2016

| Nombre | Celda |
|---|---|
| `hexagon` | Celda triangular de Chen et al. 2021 en teselado hexagonal |
| `hexagonPrint` | La misma, recortada con bordes rectos y cortes de 0.8 mm para imprimir |
| `squareThick` | Cuadrada de Rafsanjani & Pasini 2016 con t/l = 0.05 |
| `squareThin` | Cuadrada de Rafsanjani & Pasini 2016 con t/l = 0.025 |
| `curvedBeam` | Celda de viga curva importada de GiD (antes RVE) |
| `frame` | Marco hueco con punzón |

El número final es el número de celdas por lado (1 o 3).

Probetas finitas 2D con el módulo hiperelástico de Swan: neo-Hookeano con `μ = 1`
y `λ* ≈ 1.64`, elementos P1 y `ε` impuesta con pasos de 0.005. Los resultados vigentes
son los de contacto IPC (sufijo `_barrier`), con Newton amortiguado, búsqueda de línea,
curvatura negativa y predictor tangente. Los resultados sin contacto de la tabla antigua
se obtuvieron con Newton puro y en la cuadrada fina siguen una rama simétrica inestable.

## Qué imprimir en la Bambu A1

**`results/stl/squareThin_3x3_print.stl`**: 240 × 240 × 3 mm, corte 0.5 mm y bisagra
1 mm. Con IPC sigue siendo biestable en 3×3 (máximo de E en ε = 0.215, segundo mínimo en
0.305, E2/Emáx = 0.78), pero **la celda sola no lo es** (fuerza mínima ≈ 0.03 > 0). El
1×1 sirve para probar la impresión, no la biestabilidad.

Es el único caso con desplazamiento impuesto en un lado, sin solapes y biestable:
`u_y` impuesto arriba, `u_y = 0` abajo y `u_x` libre en los dos lados, salvo un nodo. Sale
de `squareThin3_UniaxialRollers_r1`, con la geometría a la mitad de tamaño; el modelo 2D no
depende de la escala.

- Máximo de energía en ε = 0.195 y segundo mínimo en ε = 0.305.
- En ε = 0.31 mide 157.2 × 157.2 mm desde 120 × 120: se expande igual en las dos direcciones.
- Estiramiento máximo en las bisagras: 1.24 en ε = 0.19, 1.32 en ε = 0.31, 1.49 en ε = 0.35.

Descartados para esta condición (1 celda, sin 3×3):

- Chen recortado con bordes rectos (`hexagonPrint`, corte 0.8 y bisagra 1.2 mm), con rodillos:
  solapa de ε = 0.07 a 0.205 y Newton para.
- Cuadrada con corte 0.8 y `t/l = 0.03`, con rodillos: solapa y no es biestable. Con los
  bordes empotrados en x: solapa y estira ×2.5.
- Cuadrada con corte 0.6 y `t/l = 0.025`, con rodillos: solapa de ε = 0.08 a 0.17 y Newton
  para en 0.24. Hay sensibilidad al ancho de corte: no cambiar el STL sin simularlo.

## Resultados

La fuerza es `dE/dε / μ`. "Biestable" significa que la fuerza cruza a negativo y
vuelve a cero, con las condiciones de contorno impuestas.

| Caso | Condición | Nodos | Llega a ε | Solapes | Biestable | Máx. E en ε | 2.º mín. en ε | E2/Emáx |
|---|---|---|---|---|---|---|---|---|
| Hexágono 1 celda | BiaxialVertices | 11 246 | 0.40 | no | sí | 0.185 | 0.26 | 0.92 |
| Hexágono 1 celda, malla fina | BiaxialVertices | 40 658 | 0.40 | no | sí | 0.185 | 0.265 | 0.90 |
| Hexágono 1 celda | UniaxialVertices | 11 246 | 0.40 | sí, ε = 0.01 a 0.125 | sí | 0.185 | 0.26 | 0.92 |
| Hexágono 1 celda | UniaxialRollers | 11 246 | 0.05 (Newton para) | no | no alcanzado | | | |
| **Hexágono 3×3** | **BiaxialVertices** | 96 212 | 0.35 | no | **sí** | 0.175 | 0.28 | 0.60 |
| Cuadrada t/l = 0.05, 1 celda | las tres | ~10 600 | 0.45 | no | **no** | | | |
| Cuadrada t/l = 0.05, 3×3 | UniaxialRollers | 94 138 | 0.35 | no | **no** | | | |
| Cuadrada t/l = 0.025, 1 celda | UniaxialRollers | 12 954 | 0.35 | no | sí | 0.215 | 0.29 | 0.89 |
| Cuadrada t/l = 0.025, 1 celda | BiaxialVertices | 12 954 | 0.35 | no | no (mínimo de fuerza > 0) | | | |
| Cuadrada t/l = 0.025, 1 celda, malla fina | BiaxialVertices | 46 358 | 0.35 | no | marginal (Fmín = -0.002) | 0.23 | 0.24 | 1.00 |
| **Cuadrada t/l = 0.025, 3×3** | **UniaxialRollers** | 115 086 | 0.35 | no | **sí** | 0.195 | 0.305 | 0.63 |

Referencia de Chen (Fig. 2d, celda periódica, leído a ojo): máximo en ε ≈ 0.18,
segundo mínimo en ε ≈ 0.29, E2/Emáx ≈ 0.45.

### Con contacto IPC (vigentes)

| Caso | Condición | Llega a ε | Solapes | Biestable | Máx. E en ε | 2.º mín. en ε | E2/Emáx |
|---|---|---|---|---|---|---|---|
| Hexágono 1 celda | BiaxialVertices | 0.40 | no | sí | 0.185 | 0.255 | 0.92 |
| **Hexágono 3×3** | **BiaxialVertices** | 0.35 | no | **sí** | 0.175 | 0.28 | 0.60 |
| Cuadrada t/l = 0.025, 1 celda | UniaxialRollers | 0.35 | no | **no** (Fmín ≈ 0.03) | | | |
| **Cuadrada t/l = 0.025, 3×3** | **UniaxialRollers** | 0.35 | no | **sí** | 0.215 | 0.305 | 0.78 |
| Viga curva 1 celda | ClampedFixedSides | -0.44 | no | sí, y tope de contacto desde -0.425 | | | |
| Viga curva 3×3 | ClampedFixedSides | -0.42 | no | saltos sucesivos: la fuerza cambia de signo en -0.235, -0.27, -0.37 y -0.40 | | | |

Chen sale igual que sin contacto: el contacto no se activa. En la cuadrada fina 3×3 la
probeta se abre primero con una cizalla asimétrica (ε = 0.08 a 0.26) y luego se ordena
en la configuración abierta simétrica.

Celdas sueltas añadidas con IPC (sin solapes en ningún caso):

| Caso | Condición | Llega a ε | Biestable |
|---|---|---|---|
| Hexágono 1 celda | UniaxialRollers | 0.40 (antes Newton paraba en 0.24) | marginal: la fuerza baja a -0.06 entre 0.135 y 0.155 |
| Hexágono 1 celda | ClampedEdges | 0.35 | no |
| hexagonPrint 1 celda | UniaxialRollers | 0.35 | no |
| Cuadrada t/l = 0.05, 1 celda | UniaxialRollers | 0.45 | no |
| Cuadrada t/l = 0.025, 1 celda | ClampedEdges | 0.32 de 0.35 (para) | no |
| Viga curva 1 celda | RollersFixedSides | -0.44 | sí, cambios de signo en -0.30 y -0.40 |
| Viga curva 1 celda | ClampedEdges y UniaxialRollers | -0.44 | no |

### Formulaciones de contacto probadas

Todas sobre el borde de la malla, en `src/Functionals`, y elegidas con `contact` en
`mainBistableSpecimens`:

- **IPC (`barrier`)**: barrera logarítmica sobre la distancia sin signo nodo-arista,
  energía infinita si se cruzan bordes. **Converge en todos los casos.**
- **Penalización (`penalty`)** nodo-arista con ventana suavizada: se atasca.
- **Lagrangiano aumentado (`augmentedLagrangian`)**, Uzawa sobre la penalización: las
  derivadas están verificadas, pero en `squareThin1` para en ε = 0.04 con el residuo
  oscilando en 1e-2, igual que la penalización. El problema es la geometría nodo-arista
  con signo (corte en profundidad, extremos de ventana, exclusión de referencia), no cómo
  se impone la restricción.
- **Tercer medio con level set**: pendiente, para hacerlo juntos.

## Marco con punzón

Marco de 36 × 9 mm con paredes de 2 mm (`frame1`, `ClampedBottomPunch`): base empotrada y
`u_y` impuesto hacia abajo en un tramo de 2 mm en el centro de arriba.

- **IPC con `d̂ = 0.4 mm` y continuación en κ** (`gapTolerance = 0.02`, κ se divide entre 10
  hasta que el hueco baja de la tolerancia, y se reinicia en cada paso): llega a 6.5 mm sin
  solapes y deja un hueco de milésimas de mm. Con κ fijo el hueco era 0.16 mm.

## Tercer medio con level set

- **Geometría:** el sólido se define con un level set sobre una malla de fondo de toda la
  caja. En el marco es `GeometricalFunction` `RectangleInclusion` sobre `TriangleMesh`. En las
  celdas es la distancia con signo a la polyshape sobre una malla de Gmsh de la caja,
  refinada junto al contorno. Después `LevelSet`, `CharacteristicFunction` y filtro `LUMP`
  dan ρ en P1.
- **Material** (`ThirdMediumFunctional`): neo-Hookeano con μ y κ interpolados con
  `SimpAllInterpolator` entre γ·sólido (ρ = 0) y sólido (ρ = 1), más HuHu
  `k_r(1-ρ)/2 |∇∇u|²` con `u` en P2. Continuación en γ dentro de cada paso: de 1e-2 al γ
  final, dividiendo entre 10.
- **Cambios en Swan:** `Hess(u)` y segundas derivadas de `TriangleQuadratic` (verificado
  con un campo cuadrático en malla distorsionada, error 1e-13). Orden de `u` configurable en
  `HyperelasticityComputer`, `DirichletCondition` y `BoundaryConditions`. El trabajo externo
  solo se calcula si hay cargas.
- **Resultado en el marco:** con γ final = 1e-6 converge hasta 4.0-4.2 mm y luego Newton
  para, porque el tercer medio bajo el punzón llega a J ≈ 0.01. No depende de `k_r`
  (1e-4, 1e-3 y 1e-2 dan lo mismo) ni del paso (0.25 o 0.05). El filtro LUMP deja la fila
  de la interfaz con ρ = 0.5, que con SimpAll es casi tan blanda como el tercer medio.
- **Rigidez final del tercer medio:** con γ final = 1e-4 llega a 5.3 mm. **Con γ final = 1e-3
  llega a 6.5 mm** y transmite el contacto (fuerza final 3.0 frente a 5.1 con IPC; el
  tercer medio y la interfaz blanda ablandan la respuesta).
- **Viga curva 1 celda con tercer medio** (`ClampedWalledSides`, γ final = 1e-6, k_r = 1e-2,
  paso -0.01): llega a ε = -0.40, que es el inicio del contacto, con 5-6 iteraciones por paso.

  | | Tercer medio | IPC |
  |---|---|---|
  | Fuerza máxima | 0.87 en 0.07 | 1.19 en 0.07 |
  | Cruce por cero | ≈ 0.295 | ≈ 0.30 |
  | Valle | -0.16 en 0.35 | -0.24 en 0.36 |
  | Contacto | ≈ 0.395 | ≈ 0.41-0.42 |

  El contacto aparece en las esquinas laterales, igual que en IPC. Después de -0.40 el medio
  llega a J = 0.007 y Newton para. k_r no cambia la fuerza; γ = 1e-5 y 1e-4 frenan el salto
  porque el medio de las cavidades se comprime antes del contacto. Lo que lo hizo funcionar:
  - **Paredes laterales** (`ClampedWalledSides`: `u_x = 0` en x = ±xmax de la caja). Sin
    ellas el medio se escapa por las aberturas laterales y para en -0.22.
  - **Amortiguamiento proximal diagonal** en `DisplacementUpdater`: `K + η·diag(|diag K|)`.
  - **Continuación en γ solo en el primer paso.** Reiniciar a γ = 1e-2 en cada paso empujaba
    la viga hacia atrás en pleno salto y paraba en -0.33.

  Pendiente: la fuerza máxima sale un 27 % por debajo de IPC (posible adelgazamiento de las
  vigas por el filtro LUMP, sin comprobar). Esta tabla salió de una copia de pruebas; falta
  repetirla con el código de la rama.
- **Mallas de fondo de las celdas:** RVE 29k elementos, hexagonPrint 128k, cuadrada 221k. Las
  geometrías de tamaño paper (cortes de 0.125 y 0.25 mm) no son abordables así.

## Solver

- **Subdivisión adaptativa del paso de carga** (`maxBisections`, 4 en el experimento): si
  Newton no converge o el punto de partida es inválido, el incremento se parte por la
  mitad y se reintenta.
- **Factorización protegida:** si la matriz o el residuo tienen valores no finitos, el
  paso falla en lugar de buscar un desplazamiento infinito. Antes, hexagon3 con rodillos se
  colgaba pidiendo una matriz llena de 550 GB.

## Condiciones de contorno

- **BiaxialVertices:** `u = ε(x - xc)` en parches pequeños en los vértices de la red del
  borde. Esos vértices están sobre líneas de espejo, así que en la expansión simétrica
  solo se mueven de forma afín y la condición acompaña al mecanismo.
- **UniaxialRollers:** `u_y` impuesto en los lados superior e inferior, `u_x` libre y un
  solo nodo con `u_x` bloqueado.
- **UniaxialVertices:** igual, pero solo en los vértices de la red de esos lados.
- **BiaxialCorners:** solo en las esquinas. Falló en ε = 0.04 antes del predictor y no se
  ha repetido.

Antes se probó empotrar el lado inferior, entero o solo la franja izquierda, y había
solapes grandes. Esos casos están en `results/fullBottomClamp` y `results/leftBottomClamp`.

## Limitaciones

- La biestabilidad se mide con el contorno impuesto. **No se ha simulado soltar la
  probeta ni la vuelta.**
- Las posiciones de Chen salen cerca del paper, pero la barrera es menor (E2/Emáx de 0.60
  en 3×3 frente a ≈ 0.45). Son probetas finitas, no periódicas, con P1 y `λ*` en lugar de
  tensión plana exacta.
- La comprobación de solapes cuenta bordes del contorno cruzados.
- El mínimo y máximo de `J` y de `tr C` se evalúan en el centroide de cada elemento en
  lugar de con `project`, que se quedaba sin memoria en las mallas grandes. En P1 es lo
  mismo (comprobado a 1e-15).
- Falta estudio de malla en 3×3 y en la cuadrada fina con rodillos.
- `squareThin3` con BiaxialVertices terminó los 61 pasos y se cortó sin error antes de
  guardar. No se ha repetido, porque la 3×3 con rodillos ya da esta celda.

## STL (`results/stl`, en mm)

Todos cerrados y con el volumen esperado. `stl_overview.png` muestra los de tamaño paper.

| Caso | `paper` | `print` | Escala | Corte / bisagra en `print` | Tamaño 3×3 en `print` |
|---|---|---|---|---|---|
| `hexagon` | 2.3 mm | 3 mm | ×4 | 0.5 / 1 mm | 200 × 242.5 mm |
| `squareThin` | 2.3 mm | 3 mm | ×2 | 0.5 / 1 mm | 240 × 240 mm |
| `squareThick` | 2.3 mm | 3 mm | ×2 | 0.5 / 2 mm | 240 × 240 mm |

- **Para imprimir, usar `print`.** Los cortes de `paper` (0.125 y 0.25 mm) están por
  debajo de lo que resuelve una impresora FDM normal.
- **`squareThick` (t/l = 0.05) no ha salido biestable en simulación.** Si se busca
  biestabilidad, imprimir `squareThin`.
- Las esquinas de los cortes son vivas. El grosor de `print` y el de Rafsanjani
  (que no lo da el paper) son decisiones mías. La escala no se ha probado en impresión.

## Supuestos

- `l0 = 10 mm` en Chen (de su Fig. 7).
- `s = 0.25 mm` y bisagra en la punta igual a `t` en Rafsanjani.

## Cambio en Swan

Con opciones que por defecto no cambian nada:

- `HyperelasticityComputer.m`: copia de `uFun` por paso, `stopWhenNotConverged`,
  predictor tangente, conteo de solapes y actualización de multiplicadores de contacto.
- `DisplacementUpdater.m`: con `lineSearch`, Newton con desplazamiento del Hessiano
  (Cholesky con reordenación), dirección de curvatura negativa, amortiguamiento
  adaptativo y convergencia por residuo. Sin `lineSearch`, el Newton de siempre.
- `ElasticityFunctional.m`: término de contacto opcional.
- Nuevos: `BoundaryContactFunctional`, `SelfContactFunctional`,
  `BarrierContactFunctional`, `AugmentedLagrangianContactFunctional`.

`HyperelasticityTests` pasa los dos casos. Sigue pendiente que `DisplacementUpdater`
calcula la reacción como `-K*u`, que no es correcta en no lineal; aquí la fuerza sale de
la fuerza interna.

## Cómo repetirlo

Desde la raíz de Swan, con `src` y esta carpeta en el path:

```matlab
runs = {'hexagon3','BiaxialVertices',0.35}; refinement = 1;
mainBistableSpecimens
mainBistableStl
```

```bash
python3 -I plotBistableCurves.py results/nombre/nombre
python3 -I renderVtuGif.py /ruta/absoluta/results/nombre/nombre 0.005
```

Cada caso guarda en `results/nombre/`: `.mat`, un `.vtu` por paso, `_setup.png`,
`_curves.png`, `_open.gif` y `_open_frames.png`.
