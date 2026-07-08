#!/bin/bash
# ==============================================================================
# SCRIPT DE BENCHMARKING PARA METAGENÓMICA ONT
# Descripción: Ejecución paralela y evaluación de herramientas alternativas.
# Uso: ./benchmark_comparativas.sh -i <lecturas_crudas.fastq.gz> -t <hilos>
# ==============================================================================

INPUT_FASTQ=$2
THREADS=$4

# Activación del entorno base para Conda
source $(conda info --base)/etc/profile.d/conda.sh

mkdir -p BENCHMARKS/{Fase1_Filtrado,Fase2_Ensamblaje,Fase3_Pulido,Fase4_Binning,Fase5_Anotacion}/{resultados,evaluacion}

# ==============================================================================
# FASE 1: BENCHMARK DE CONTROL DE CALIDAD Y FILTRADO
# ==============================================================================
echo "[FASE 1] Iniciando comparativa: Filtlong vs Chopper..."
conda activate ont_prueba1

# Evaluación del control (Datos crudos)
NanoStat --fastq $INPUT_FASTQ --threads $THREADS -n BENCHMARKS/Fase1_Filtrado/evaluacion/stats_crudo.txt

# Ejecución de Filtlong
filtlong --min_length 1000 --min_mean_q 9 $INPUT_FASTQ > BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq
NanoStat --fastq BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq --threads $THREADS -n BENCHMARKS/Fase1_Filtrado/evaluacion/stats_filtlong.txt

# Ejecución de Chopper
conda activate ont_prueba2
gunzip -c $INPUT_FASTQ | chopper -q 9 -l 1000 > BENCHMARKS/Fase1_Filtrado/resultados/chopper_clean.fastq
NanoStat --fastq BENCHMARKS/Fase1_Filtrado/resultados/chopper_clean.fastq --threads $THREADS -n BENCHMARKS/Fase1_Filtrado/evaluacion/stats_chopper.txt

# ==============================================================================
# FASE 2: BENCHMARK DE ENSAMBLAJE METAGENÓMICO
# ==============================================================================
echo "[FASE 2] Iniciando comparativa: Flye vs Raven..."
conda activate ont_prueba1

# Ejecución de Flye
flye --nano-raw BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq --meta --out-dir BENCHMARKS/Fase2_Ensamblaje/resultados/flye_out --threads $THREADS

# Ejecución de Raven
conda activate ont_prueba2
raven --threads $THREADS BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq > BENCHMARKS/Fase2_Ensamblaje/resultados/raven_assembly.fasta

# Evaluación comparativa con MetaQUAST
conda activate juez_quast
metaquast.py BENCHMARKS/Fase2_Ensamblaje/resultados/flye_out/assembly.fasta BENCHMARKS/Fase2_Ensamblaje/resultados/raven_assembly.fasta -o BENCHMARKS/Fase2_Ensamblaje/evaluacion/metaquast_out --max-ref-number 0 --threads $THREADS

# ==============================================================================
# FASE 3: BENCHMARK DE PULIDO ESTRUCTURAL
# ==============================================================================
echo "[FASE 3] Iniciando comparativa de iteraciones de pulido: Racon 1x vs Racon 2x..."
conda activate ont_prueba2

# Racon 1 Iteración
minimap2 -x map-ont -t $THREADS BENCHMARKS/Fase2_Ensamblaje/resultados/raven_assembly.fasta BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq > BENCHMARKS/Fase3_Pulido/resultados/map1.paf
racon -t $THREADS BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq BENCHMARKS/Fase3_Pulido/resultados/map1.paf BENCHMARKS/Fase2_Ensamblaje/resultados/raven_assembly.fasta > BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta

# Racon 2 Iteraciones
minimap2 -x map-ont -t $THREADS BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq > BENCHMARKS/Fase3_Pulido/resultados/map2.paf
racon -t $THREADS BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq BENCHMARKS/Fase3_Pulido/resultados/map2.paf BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta > BENCHMARKS/Fase3_Pulido/resultados/raven_racon2.fasta

# Evaluación del impacto con BUSCO
echo "[EVALUACIÓN 3] Ejecutando BUSCO para cuantificar la integridad de los marcos de lectura..."
busco -i BENCHMARKS/Fase2_Ensamblaje/resultados/raven_assembly.fasta -l bacteria_odb10 -o busco_crudo -m genome -c $THREADS --out_path BENCHMARKS/Fase3_Pulido/evaluacion/
busco -i BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta -l bacteria_odb10 -o busco_racon1x -m genome -c $THREADS --out_path BENCHMARKS/Fase3_Pulido/evaluacion/
busco -i BENCHMARKS/Fase3_Pulido/resultados/raven_racon2.fasta -l bacteria_odb10 -o busco_racon2x -m genome -c $THREADS --out_path BENCHMARKS/Fase3_Pulido/evaluacion/

# ==============================================================================
# FASE 4: BENCHMARK DE RECONSTRUCCIÓN DE GENOMAS (BINNING)
# ==============================================================================
echo "[FASE 4] Preparando perfil de cobertura para MetaBAT2 y CONCOCT..."
conda activate ont_prueba1

mkdir -p BENCHMARKS/Fase4_Binning/{alineamiento,metabat2,concoct_out,concoct_bins,dastool}
minimap2 -ax map-ont -t $THREADS BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta BENCHMARKS/Fase1_Filtrado/resultados/filtlong_clean.fastq | samtools sort -@ $THREADS -o BENCHMARKS/Fase4_Binning/alineamiento/cobertura.bam
samtools index BENCHMARKS/Fase4_Binning/alineamiento/cobertura.bam

echo "[FASE 4] Ejecución de MetaBAT2..."
jgi_summarize_bam_contig_depths --outputDepth BENCHMARKS/Fase4_Binning/metabat2/raven_depth.txt --percentIdentity 85 --minContigLength 1500 --minContigDepth 1.0 --referenceFasta BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta BENCHMARKS/Fase4_Binning/alineamiento/cobertura.bam
metabat2 -i BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta -a BENCHMARKS/Fase4_Binning/metabat2/raven_depth.txt -o BENCHMARKS/Fase4_Binning/metabat2/bin -m 1500

echo "[FASE 4] Ejecución de CONCOCT..."
conda activate ont_prueba2
cut_up_fasta.py BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta -c 10000 -o 0 --merge_last -b BENCHMARKS/Fase4_Binning/concoct_out/contigs_10k.bed > BENCHMARKS/Fase4_Binning/concoct_out/contigs_10k.fa
concoct_coverage_table.py BENCHMARKS/Fase4_Binning/concoct_out/contigs_10k.bed BENCHMARKS/Fase4_Binning/alineamiento/cobertura.bam > BENCHMARKS/Fase4_Binning/concoct_out/coverage_table.tsv
concoct --composition_file BENCHMARKS/Fase4_Binning/concoct_out/contigs_10k.fa --coverage_file BENCHMARKS/Fase4_Binning/concoct_out/coverage_table.tsv -b BENCHMARKS/Fase4_Binning/concoct_out/ -t $THREADS
merge_cutup_clustering.py BENCHMARKS/Fase4_Binning/concoct_out/clustering_gt1000.csv > BENCHMARKS/Fase4_Binning/concoct_out/clustering_merged.csv
extract_fasta_bins.py BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta BENCHMARKS/Fase4_Binning/concoct_out/clustering_merged.csv --output_path BENCHMARKS/Fase4_Binning/concoct_bins/

echo "[FASE 4] Refinamiento de consenso con DAS_Tool..."
conda activate ont_prueba1
Fasta_to_Contig2Bin.sh -i BENCHMARKS/Fase4_Binning/metabat2 -e fa > BENCHMARKS/Fase4_Binning/metabat2.tsv
Fasta_to_Contig2Bin.sh -i BENCHMARKS/Fase4_Binning/concoct_bins -e fa > BENCHMARKS/Fase4_Binning/concoct.tsv
DAS_Tool -i BENCHMARKS/Fase4_Binning/metabat2.tsv,BENCHMARKS/Fase4_Binning/concoct.tsv -l metabat,concoct -c BENCHMARKS/Fase3_Pulido/resultados/raven_racon1.fasta -o BENCHMARKS/Fase4_Binning/dastool -t $THREADS --write_bins

# Evaluación de MAGs
echo "[EVALUACIÓN 4] Ejecutando CheckM2 para control de calidad genómica..."
conda activate juez_checkm2
checkm2 predict --threads $THREADS --input BENCHMARKS/Fase4_Binning/metabat2/ -x fa --output-dir BENCHMARKS/Fase4_Binning/evaluacion/checkm2_metabat
checkm2 predict --threads $THREADS --input BENCHMARKS/Fase4_Binning/concoct_bins/ -x fa --output-dir BENCHMARKS/Fase4_Binning/evaluacion/checkm2_concoct
checkm2 predict --threads $THREADS --input BENCHMARKS/Fase4_Binning/dastool_DASTool_bins/ -x fa --output-dir BENCHMARKS/Fase4_Binning/evaluacion/checkm2_dastool


# ==============================================================================
# FASE 5: BENCHMARK DE ANOTACIÓN TAXONÓMICA Y FUNCIONAL
# ==============================================================================
echo "[FASE 5] Iniciando comparativas de anotación taxonómica y funcional..."
mkdir -p BENCHMARKS/Fase5_Anotacion/{gtdbtk,kraken2,bakta,prokka}

# Para el benchmark dinámico, iteramos sobre los MAGs resultantes de DAS_Tool
for mag in BENCHMARKS/Fase4_Binning/dastool_DASTool_bins/*.fa; do
    mag_name=$(basename "$mag" .fa)
    
    # 1. Taxonomía (GTDB-Tk vs Kraken2)
    echo ">> Clasificando $mag_name ..."
    conda activate ont_prueba1
    # Nota: GTDB-Tk procesa un directorio completo, se ejecutará fuera del bucle para mayor eficiencia.
    
    conda activate anotacion_env
    kraken2 --db ~/Doctorado/bases_datos/kraken2_db --threads $THREADS --use-names --report BENCHMARKS/Fase5_Anotacion/kraken2/reporte_${mag_name}.txt $mag > BENCHMARKS/Fase5_Anotacion/kraken2/salida_${mag_name}.kraken
    
    # 2. Función (Bakta vs Prokka)
    echo ">> Anotando funcionalmente $mag_name ..."
    bakta --db ~/Doctorado/bases_datos/bakta_db_light/db-light --threads $THREADS --skip-sorf --output BENCHMARKS/Fase5_Anotacion/bakta/${mag_name} --prefix ${mag_name} --force $mag
    prokka --outdir BENCHMARKS/Fase5_Anotacion/prokka/${mag_name} --prefix ${mag_name} --cpus $THREADS --force $mag
done

# Ejecución batch de GTDB-Tk
conda activate ont_prueba1
gtdbtk classify_wf --genome_dir BENCHMARKS/Fase4_Binning/dastool_DASTool_bins/ --out_dir BENCHMARKS/Fase5_Anotacion/gtdbtk/ --extension fa --cpus $THREADS

echo "[COMPLETADO] Todas las pruebas empíricas han concluido. Revise los subdirectorios de evaluación para consultar las métricas completas."
