# Bioinformatics-Msc-s
# Long-Read Metagenomic Pipeline for Anaerobic Bioreactors

## Descripción
Este repositorio contiene un flujo de trabajo (*pipeline*) bioinformático modular y automatizado, diseñado específicamente para el procesamiento de datos metagenómicos basados en lecturas largas procedentes de Oxford Nanopore Technologies (ONT). 

El *pipeline* abarca desde el preprocesamiento de lecturas crudas hasta la recuperación de Genomas Ensamblados a partir de Metagenomas (MAGs) de alta pureza y su posterior caracterización funcional y taxonómica. El flujo de trabajo ha sido optimizado para maximizar la contigüidad genómica y minimizar la formación de quimeras y errores estructurales derivados de la tecnología de secuenciación.

## Estructura del Flujo de Trabajo
El pipeline óptimo integrado en este repositorio ejecuta de forma secuencial las siguientes fases:
1. **Control de Calidad y Filtrado:** Retención de lecturas largas (>1000 pb) y filtrado de calidad (Q9) utilizando `Filtlong`.
2. **Ensamblaje Metagenómico *De Novo*:** Reconstrucción genómica maximizando la contigüidad (N50) mediante el algoritmo basado en OLC de `Raven`.
3. **Pulido Estructural (*Polishing*):** Corrección de *indels* mediante una única iteración de mapeo y consenso utilizando `Minimap2` y `Racon`, evitando el sobre-pulido para preservar los marcos de lectura abiertos (validados vía `BUSCO`).
4. **Binning y Refinamiento:** Agrupación probabilística mediante la combinación de `MetaBAT2` y `CONCOCT`. Los perfiles generados son integrados y depurados utilizando `DAS_Tool` para la extracción de MAGs consenso no redundantes.
5. **Evaluación de Calidad:** Cuantificación de la completitud y la contaminación genómica con algoritmos de Machine Learning (`CheckM2`).
6. **Asignación Taxonómica y Funcional:**
   - **Taxonomía:** Identificación de linajes consensuados empleando alineamiento de marcadores universales (`GTDB-Tk`) y clasificación rápida por k-mers (`Kraken2`).
   - **Función:** Anotación metabólica estructural *alignment-free* mediante `Bakta` (ejecutado sin búsqueda de sORFs para optimización de memoria), detectando CDS, ncRNAs y sistemas CRISPR.

## Dependencias y Requisitos
Se recomienda la gestión de dependencias a través de entornos virtuales `conda` para evitar conflictos. Las principales herramientas requeridas son:
* Filtlong, NanoStat
* Raven, MetaQUAST
* Minimap2, Racon, BUSCO
* Samtools, MetaBAT2, CONCOCT, DAS_Tool
* CheckM2, GTDB-Tk, Kraken2, Bakta

## Uso
El script automatizado `pipeline_optimo.sh` ejecuta el flujo completo de principio a fin, asumiendo la disposición de las bases de datos de referencia configuradas localmente.

```bash
chmod +x pipeline_optimo.sh
./pipeline_optimo.sh -i data/raw/secuencias.fastq.gz -o resultados_pipeline -t 16
