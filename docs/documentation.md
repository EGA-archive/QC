# BAM/CRAM QC Pipelines Documentation

## Overview
The European Genome-phenome Archive (EGA) currently stores nearly 3 million BAM and CRAM files, a number that continues to grow thanks to the tireless efforts of the scientific community. To enhance the quality reports we generate for each of these files, we've developed two custom pipelines that integrate various bioinformatics tools for comprehensive quality evaluation.

- `BAM_pipeline_2.py`: Computes a wide array of quality metrics.
- `BAM_finalize_2.py`: Post-processes the outputs and adapts them for visualization in MultiQC.

This documentation provides an in-depth overview of these tools, including the modules developed and the modifications made to integrate with MultiQC.

---

## CRAM to BAM Conversion
Working directly with CRAM files is often problematic due to inconsistent support across tools. We therefore convert CRAM files to BAM using `samtools view -b`, which requires the MD5 checksum of the reference genome unless the reference sequence is embedded in the CRAM file.

---

## Tools Used

### Qualimap
Generates comprehensive BAM statistics such as:
- Mapping quality
- Coverage
- Nucleotide composition
- Homopolymer presence

It also produces a summary file (`genome_results.txt`) with global metrics.

### RSeQC
Uses a BED annotation file to determine where reads fall (exons, introns, UTRs) and creates informative plots.

### Samtools (depth -a)
Provides per-base coverage in three columns: chromosome, position, and coverage.

### Picard (CollectAlignmentSummaryMetrics)
Generates alignment quality metrics. We extract `PCT_Chimeras` for inclusion in our summary report.

---

## Custom Functions in the Pipeline

### `parse_bam_header()`
Extracts metadata from the BAM header:
- Format version
- Sort order
- Platform
- Sample name
- Tools used

## `calculate_stat_from_counter()`

Esta es una función auxiliar que se usa dentro de `calculate_coverage_stat_streaming()` para calcular una estadística (por ahora, la mediana) a partir de un conjunto de valores de cobertura.

La entrada es un `Counter`, donde cada clave representa un valor de cobertura (por ejemplo, 0, 1, 2...), y el valor asociado indica cuántas veces aparece en una ventana del genoma.  
Esta función transforma ese diccionario en una lista de valores repetidos según su frecuencia, y calcula la mediana. Si no hay valores, devuelve `0.0`.

---

## `calculate_coverage_stat_streaming()`

Esta función divide el genoma en intervalos de un tamaño fijo (por ejemplo, 3 millones de pares de bases).  
Para cada intervalo o ventana, calcula un único valor de cobertura representativo, que puede ser la **media** o la **mediana** de las coberturas observadas en ese tramo.

### 📌 ¿Cómo lo hace?

La idea principal es **resumir la cobertura del BAM** en bloques de tamaño fijo, llamados **ventanas**, a lo largo de todo el genoma. En este caso, cada ventana cubre 3.000.000 pares de bases.

---

### 🔹 Paso a paso

#### 1. Lectura línea por línea con `samtools depth`

Se lanza el comando:

```
samtools depth -a archivo.bam
```

Este comando devuelve, para cada posición del genoma, tres columnas:

```
cromosoma   posición   cobertura
chr1        0          12  
chr1        1          13  
chr1        2          15  
...
```

La cobertura representa el número de lecturas que cubren esa posición.

---

#### 2. Asignar cada posición a una ventana

Cada posición se traduce en una **posición global** (sumando offsets entre cromosomas), y a partir de esa posición se determina a qué ventana pertenece.

```python
window_start = (global_pos // window_size) * window_size
```

**Ejemplo:**  
La posición global `4.200.000` pertenece a la ventana `[3.000.000 – 5.999.999]`, que empieza en `3.000.000`.

---

#### 3. Acumular la frecuencia de cada valor de cobertura en cada ventana

Se construye un diccionario `coverage_by_window` donde:

- La **clave** es el inicio de la ventana (`window_start`)
- El **valor** es un `Counter`, que almacena cuántas veces aparece cada nivel de cobertura dentro de esa ventana.

**Ejemplo simplificado:**

```
posición   cobertura
3.000.005  0
3.000.010  0
3.000.050  1
3.000.300  1
3.000.700  2
```

Entonces:

```python
coverage_by_window[3000000] = Counter({0: 2, 1: 2, 2: 1})
```

---

### 🔁 ¿Y cómo se resume todo eso?

Una vez tienes todos los contadores (`Counter`) por ventana, necesitas obtener un único valor representativo por ventana: la **media** o la **mediana** de cobertura.

Ahí entra la función `calculate_stat_from_counter(counter)`.

```python
# Entrada:
Counter({0: 2, 1: 2, 2: 1})

# Convertido a lista:
[0, 0, 1, 1, 2]

# Mediana:
1
```

Este valor se asocia a la ventana en el resultado final.

---

## 🔢 ¿Qué es una posición global?

Como la intención es representar la cobertura a lo largo de **todo el genoma como si fuera una línea recta**, cada posición necesita tener un valor único que indique su lugar en ese “genoma concatenado”.

Para lograr esto, se construye un diccionario de **offsets acumulados**, donde:

- `chr1` empieza en la posición `0`
- `chr2` empieza justo después del final de `chr1` (por ejemplo, en la posición `249.250.621`)
- `chr3` empieza después de `chr2`, y así sucesivamente

De este modo, si una ventana empieza en la posición global `3.000.000`, sabremos si corresponde a `chr1`, `chr2`, etc., dependiendo de estos offsets.  
Es una forma práctica de recorrer el genoma como si fuera una única secuencia.

---

## 🔀 Gestión de nombres de cromosomas

Una parte importante del script está dedicada a gestionar los diferentes formatos en que aparecen los nombres de los cromosomas en los archivos BAM (por ejemplo, `1` vs `chr1`, `MT` vs `chrM`, etc.).  

Como no hay un estándar único (NCBI, UCSC, Ensembl usan nomenclaturas distintas), se usa un archivo llamado `chrnames.json`, que contiene un diccionario de equivalencias.

Si el BAM tiene nombres desconocidos, el script simplemente usa los nombres del header y los incluye tal cual en el JSON final.

---

## ⚠️ Tratamiento especial de `chrM`

El cromosoma mitocondrial (`chrM`) es muy pequeño (por ejemplo, 16.569 pb), por lo que **no se agrupa por ventanas**.  

En lugar de eso, se calcula un **valor único de cobertura** (media o mediana) para él y se coloca al final del JSON con una **posición artificial**.

---

## 🧾 Resultado final

La salida es un archivo `.json` con el siguiente formato:

```json
[
  ["chr1", 0, 24.2],
  ["chr1", 3000000, 23.8],
  ["chr2", 249250621, 27.5],
  ...
  ["chrM", 900000000, 98.1]
]
```

Cada línea representa:

- El **cromosoma** al que pertenece la ventana  
- La **posición global** en el genoma (calculada con offsets)  
- El **valor de cobertura** en esa ventana (media o mediana)


