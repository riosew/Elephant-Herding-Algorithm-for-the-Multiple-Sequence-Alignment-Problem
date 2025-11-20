# EHO – Elephant Herding Optimization for Multiple Sequence Alignment (MSA)
## English

### Overview

This repository contains an implementation of the **Elephant Herding Optimization (EHO)** algorithm applied to the problem of **Multiple Sequence Alignment (MSA)** of protein or DNA sequences.

The approach is inspired by the **herding behavior and clan organization of elephants**, where each elephant represents a candidate alignment solution, and the population evolves towards better solutions using bio-inspired operators.

The main code is in `EHO.py`.

### Algorithm features

#### Solution representation (elephants)

Each **Elephant** represents a multiple alignment of the input sequences.

Internally, an elephant stores:

- A set of aligned sequences (with `-` gaps when needed).
- A **fitness** value that measures the quality of the alignment.

Sequences are read from a **FASTA** file and are padded to the same length using the `-` character when necessary.

#### Clan behavior and herding

The algorithm follows the **elephant herding** metaphor:

- The population is organized into one or more **clans**.
- Within each clan, some elephants move towards better positions (better alignments), using fitter elephants as references.
- Others are replaced or relocated, simulating **dispersion** and **exploration** of the search space.
- This process is repeated for a fixed number of iterations, progressively refining the alignment quality.

Overall, the algorithm combines:

- **Exploration**: generating new alignment configurations.
- **Exploitation**: improving already promising solutions within the clans.

### Scoring with BLOSUM matrix and gap penalty

Alignment quality is evaluated using a **BLOSUM** substitution matrix (e.g., BLOSUM62), which assigns a score to each possible pair of aligned amino acids.

For a given alignment, the fitness function:

- Sums BLOSUM scores over the alignment columns across all sequences.
- Applies a **gap penalty** (e.g., `GAP_PENALTY = -2`) each time a `-` gap appears in an aligned position.

Thus, alignments that produce reasonable substitutions according to BLOSUM and minimize gaps tend to achieve higher scores.

The code uses the `blosum` Python package or, if not available, a built-in BLOSUM matrix defined in the script.

### Usage

1. Clone this repository or download the files.
2. Make sure you have the required dependencies installed (e.g., `numpy` and `blosum`).
3. Prepare your FASTA file containing the sequences to be aligned.
4. In `EHO.py`, replace the string `yourFastaFile.fasta` with the path to your FASTA file.

Example:

```python
# Before
sequences = read_fasta("yourFastaFile.fasta")

# After
sequences = read_fasta("mySequences.fasta")
```

5. Run `EHO.py` with Python.

```bash
python EHO.py
```

### Input format (FASTA)

The FASTA file must follow the standard format:

```text
>seq1
MKT...
>seq2
MRT...
>seq3
MKL...
```

Each sequence starts with a header line beginning with `>`, and the following lines up to the next `>` (or end of file) correspond to the sequence.

The algorithm:

- Reads the sequences.
- Removes line breaks.
- Pads with `-` up to the maximum length when needed.

### Output

The typical output includes:

- The **best alignment** found for all sequences.
- The associated **fitness score**.
- Additional information such as the number of fitness evaluations (`NFE`) and approximate runtime, depending on the code version.

You can adapt the final part of `EHO.py` in order to:

- Save the alignment to a file.
- Print it in a specific format.
- Integrate it with other bioinformatics tools.

### Citation and bibliographic reference

If you use this code in academic work or publications, please cite the article associated with this algorithm.

**Reference (to be completed manually by the repository author):**

> [Full scientific article reference associated with the EHO algorithm for MSA will go here]

---

## Licencia / License

Este proyecto puede distribuirse bajo la licencia que el autor elija (por ejemplo, MIT, GPL, etc.).

This project can be distributed under any license chosen by the author (e.g., MIT, GPL, etc.).

## Español

### Descripción general

Este repositorio contiene una implementación del algoritmo **Elephant Herding Optimization (EHO)** aplicado al problema de **alineamiento múltiple de secuencias (MSA)** de proteínas o ADN.

El enfoque está inspirado en el **comportamiento de pastoreo y organización en clanes de elefantes**, donde cada elefante representa una posible solución de alineamiento, y el grupo evoluciona hacia mejores soluciones mediante operadores bioinspirados.

El código principal se encuentra en el archivo `EHO.py`.

### Características del algoritmo

#### Representación de soluciones (elefantes)

Cada **Elephant** representa un alineamiento múltiple de las secuencias de entrada.

Internamente, un elefante almacena:

- Un conjunto de secuencias alineadas (con huecos `-` cuando es necesario).
- Un valor de **fitness**, que mide la calidad del alineamiento.

Las secuencias se leen desde un archivo en formato **FASTA** y se normalizan a la misma longitud mediante un relleno con el carácter `-` cuando es necesario.

#### Comportamiento de clanes y pastoreo

El algoritmo sigue la metáfora del **pastoreo de elefantes**:

- La población se organiza en uno o varios **clanes**.
- Dentro de cada clan, algunos elefantes se mueven hacia mejores posiciones (mejores alineamientos), tomando como referencia elefantes más aptos.
- Otros son reemplazados o reubicados, simulando la **dispersión** y la **exploración del espacio de búsqueda**.
- Este proceso se repite durante un número determinado de iteraciones, refinando progresivamente la calidad del alineamiento.

En términos generales, el algoritmo combina:

- **Exploración**: generación de nuevas configuraciones de alineamiento.
- **Explotación**: mejora de soluciones ya prometedoras dentro de los clanes.

### Evaluación con matriz BLOSUM y penalización por huecos

La calidad de cada alineamiento se evalúa usando una matriz de sustitución **BLOSUM** (por ejemplo, BLOSUM62), que asigna una puntuación a cada posible par de aminoácidos alineados.

Para un alineamiento dado, la función de fitness:

- Suma las puntuaciones de BLOSUM para las columnas del alineamiento entre todas las secuencias.
- Aplica una **penalización por huecos** (gap penalty), típicamente un valor negativo como `GAP_PENALTY = -2`, cada vez que aparece un hueco `-` en una posición alineada.

De este modo, alineamientos que conservan sustituciones razonables según BLOSUM y minimizan huecos tienden a obtener puntuaciones más altas.

El código hace uso de la biblioteca `blosum` o, en su defecto, de una matriz BLOSUM codificada internamente si no se encuentra la matriz externa.

### Uso

1. Clona este repositorio o descarga los archivos.
2. Asegúrate de tener instaladas las dependencias necesarias (por ejemplo, `numpy` y `blosum`).
3. Prepara tu archivo FASTA con las secuencias que quieres alinear.
4. En el código `EHO.py`, reemplaza la cadena `yourFastaFile.fasta` por la ruta a tu archivo FASTA.

Ejemplo:

```python
# Antes
sequences = read_fasta("yourFastaFile.fasta")

# Después
sequences = read_fasta("misSecuencias.fasta")
```

5. Ejecuta `EHO.py` con Python.

```bash
python EHO.py
```

### Formato de entrada (FASTA)

El archivo FASTA debe seguir el formato estándar:

```text
>seq1
MKT...
>seq2
MRT...
>seq3
MKL...
```

Cada secuencia comienza con una línea que inicia con `>` y las líneas siguientes hasta la próxima cabecera `>` (o fin de archivo) corresponden a la secuencia.

El algoritmo:

- Lee las secuencias.
- Elimina saltos de línea.
- Rellena con `-` hasta la longitud máxima del conjunto cuando es necesario.

### Salida

La salida típica incluye:

- El **mejor alineamiento** encontrado para todas las secuencias.
- La **puntuación de fitness** asociada.
- Información adicional como el número de evaluaciones de la función de fitness (`NFE`) y el tiempo de ejecución aproximado, dependiendo de la versión del código.

Puedes adaptar la parte final de `EHO.py` para:

- Guardar el alineamiento en un archivo.
- Imprimirlo con un formato específico.
- Integrarlo con otras herramientas de análisis bioinformático.

### Cita y referencia bibliográfica

Por favor, si utilizas este código en trabajos académicos o publicaciones, cita el artículo asociado a este algoritmo.

**Referencia (a completar manualmente por el autor del repositorio):**

> [Aquí irá la referencia completa del artículo científico asociado al algoritmo EHO para MSA]


---
