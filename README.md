# Bioinformatics-Msc-s
Repository containing the projects and code developed during my Master’s in Bioinformatics. Includes analyses, scripts, and exercises focused on computational biology and biological data processing.
# Metagenomic Benchmarks: Comparative Analysis of Long-Read Tools

## Descripción
Este repositorio complementario contiene los scripts y la documentación necesaria para reproducir empíricamente las comparativas metodológicas realizadas durante el diseño del pipeline metagenómico principal. 

El objetivo de este módulo es garantizar la total transparencia y reproducibilidad de las decisiones algorítmicas tomadas a lo largo del estudio, permitiendo a los usuarios ejecutar los mismos *benchmarks* y obtener las mismas métricas de evaluación para sus propios conjuntos de datos procedentes de Oxford Nanopore Technologies (ONT).

## Comparativas Evaluadas
El script de validación `benchmark_comparativas.sh` ejecuta los siguientes enfrentamientos algorítmicos paralelos:
1. **Filtro de Calidad:** `Filtlong` (evaluación holística) vs. `Chopper` (procesamiento lineal en flujo). Evaluado mediante estadísticos de contigüidad y calidad con `NanoStat`.
2. **Ensamblaje *De Novo*:** `Flye` (grafos de repetición) vs. `Raven` (Overlap-Layout-Consensus). Evaluado mediante métricas de continuidad (N50) y fragmentación estructural con `MetaQUAST`.
3. **Pulido Estructural:** Ensamblaje Crudo vs. `Racon` (1 iteración) vs. `Racon` (2 iteraciones). Evaluado mediante la restauración de marcos de lectura y genes ortólogos de copia única con `BUSCO`.
4. **Anotación Funcional y Taxonómica:** Contrastes de rendimiento entre `GTDB-Tk` y `Kraken2` (taxonomía), y `Bakta` frente a `Prokka` (anotación metabólica estructural).

## Requisitos Previos
Para la correcta ejecución del benchmarking, se requiere la instalación de las herramientas de evaluación asociadas:
* NanoStat, MetaQUAST, BUSCO (linaje `bacteria_odb10`).

## Instrucciones de Uso
Para iniciar la batería de pruebas comparativas, ejecute el script proporcionado indicando el archivo de lecturas crudas y los hilos de procesamiento disponibles:

```bash
chmod +x benchmark_comparativas.sh
./benchmark_comparativas.sh -i data/01_raw/ERR5363646.fastq.gz -t 16
