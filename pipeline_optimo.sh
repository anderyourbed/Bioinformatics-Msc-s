### 2. Archivo `pipeline_optimo.sh` (El Script)
```bash
#!/bin/bash
# ==============================================================================
# SCRIPT ÓPTIMO DE METAGENÓMICA (ONT LONG-READS)
# Descripción: Flujo automatizado de filtrado, ensamblaje, pulido, binning y anotación.
# Uso: ./pipeline_optimo.sh <input_fastq.gz> <output_dir> <threads>
# ==============================================================================

# Argumentos de entrada
INPUT_FASTQ=$1
OUT_DIR=$2
THREADS=$3

if [ -z "$INPUT_FASTQ" ] || [ -z "$OUT_DIR" ] || [ -z "$THREADS" ]; then
    echo "Uso: ./pipeline_optimo.sh <input_fastq.gz> <directorio_salida> <hilos>"
    exit 1
fi

# Habilitar activación de entornos conda dentro del script
source $(conda info --base)/etc/profile.d/conda.sh

# Creación de estructura de directorios
mkdir -p $OUT_DIR/{01_qc,02_assembly,03_polishing,04_binning,05_evaluation,06_annotation}

# ------------------------------------------------------------------------------
# FASE 1: CONTROL DE CALIDAD Y FILTRADO (FILTLONG)
# ------------------------------------------------------------------------------
echo "[FASE 1] Ejecutando filtrado de lecturas (Longitud > 1000pb, Calidad > Q9)..."
conda activate ont_prueba1

filtlong --min_length 1000 --min_mean_q 9 $INPUT_FASTQ > $OUT_DIR/01_qc/clean_reads.fastq

# ------------------------------------------------------------------------------
# FASE 2: ENSAMBLAJE METAGENÓMICO (RAVEN)
# ------------------------------------------------------------------------------
echo "[FASE 2] Ensamblando metagenoma con Raven..."
conda activate ont_prueba2

raven --threads $THREADS $OUT_DIR/01_qc/clean_reads.fastq > $OUT_DIR/02_assembly/raven_assembly.fasta

# ------------------------------------------------------------------------------
# FASE 3: PULIDO ESTRUCTURAL - 1 ITERACIÓN (MINIMAP2 + RACON)
# ------------------------------------------------------------------------------
echo "[FASE 3] Aplicando pulido estructural de consenso (Racon 1x)..."
minimap2 -x map-ont -t $THREADS $OUT_DIR/02_assembly/raven_assembly.fasta $OUT_DIR/01_qc/clean_reads.fastq > $OUT_DIR/03_polishing/map.paf

racon -t $THREADS $OUT_DIR/01_qc/clean_reads.fastq $OUT_DIR/03_polishing/map.paf $OUT_DIR/02_assembly/raven_assembly.fasta > $OUT_DIR/03_polishing/assembly_polished.fasta

# ------------------------------------------------------------------------------
# FASE 4: BINNING Y REFINAMIENTO (METABAT2 + CONCOCT + DAS_TOOL)
# ------------------------------------------------------------------------------
echo "[FASE 4] Mapeando lecturas para perfil de cobertura..."
conda activate ont_prueba1
minimap2 -ax map-ont -t $THREADS $OUT_DIR/03_polishing/assembly_polished.fasta $OUT_DIR/01_qc/clean_reads.fastq | samtools sort -@ $THREADS -o $OUT_DIR/04_binning/cobertura.bam
samtools index $OUT_DIR/04_binning/cobertura.bam

echo "[FASE 4] Ejecutando MetaBAT2..."
jgi_summarize_bam_contig_depths --outputDepth $OUT_DIR/04_binning/depth.txt --percentIdentity 85 --minContigLength 1500 --minContigDepth 1.0 --referenceFasta $OUT_DIR/03_polishing/assembly_polished.fasta $OUT_DIR/04_binning/cobertura.bam
metabat2 -i $OUT_DIR/03_polishing/assembly_polished.fasta -a $OUT_DIR/04_binning/depth.txt -o $OUT_DIR/04_binning/metabat2/bin -m 1500

echo "[FASE 4] Ejecutando CONCOCT..."
conda activate ont_prueba2
mkdir -p $OUT_DIR/04_binning/concoct
cut_up_fasta.py $OUT_DIR/03_polishing/assembly_polished.fasta -c 10000 -o 0 --merge_last -b $OUT_DIR/04_binning/concoct/contigs_10k.bed > $OUT_DIR/04_binning/concoct/contigs_10k.fa
concoct_coverage_table.py $OUT_DIR/04_binning/concoct/contigs_10k.bed $OUT_DIR/04_binning/cobertura.bam > $OUT_DIR/04_binning/concoct/coverage_table.tsv
concoct --composition_file $OUT_DIR/04_binning/concoct/contigs_10k.fa --coverage_file $OUT_DIR/04_binning/concoct/coverage_table.tsv -b $OUT_DIR/04_binning/concoct/ -t $THREADS
merge_cutup_clustering.py $OUT_DIR/04_binning/concoct/clustering_gt1000.csv > $OUT_DIR/04_binning/concoct/clustering_merged.csv
extract_fasta_bins.py $OUT_DIR/03_polishing/assembly_polished.fasta $OUT_DIR/04_binning/concoct/clustering_merged.csv --output_path $OUT_DIR/04_binning/concoct/bins/

echo "[FASE 4] Refinamiento de consenso con DAS_Tool..."
conda activate ont_prueba1
Fasta_to_Contig2Bin.sh -i $OUT_DIR/04_binning/metabat2 -e fa > $OUT_DIR/04_binning/metabat2.tsv
Fasta_to_Contig2Bin.sh -i $OUT_DIR/04_binning/concoct/bins -e fa > $OUT_DIR/04_binning/concoct.tsv
DAS_Tool -i $OUT_DIR/04_binning/metabat2.tsv,$OUT_DIR/04_binning/concoct.tsv -l metabat,concoct -c $OUT_DIR/03_polishing/assembly_polished.fasta -o $OUT_DIR/04_binning/dastool -t $THREADS --write_bins

# ------------------------------------------------------------------------------
# FASE 5: EVALUACIÓN DE PUREZA (CHECKM2)
# ------------------------------------------------------------------------------
echo "[FASE 5] Evaluando la completitud y contaminación de los MAGs consenso..."
conda activate juez_checkm2
checkm2 predict --threads $THREADS --input $OUT_DIR/04_binning/dastool_DASTool_bins/ -x fa --output-dir $OUT_DIR/05_evaluation/

# ------------------------------------------------------------------------------
# FASE 6: ANOTACIÓN TAXONÓMICA Y FUNCIONAL
# ------------------------------------------------------------------------------
# Se toma el mejor MAG generado por DAS_Tool (ejemplo: bin.3.fa, requiere automatización para seleccionarlo dinámicamente)
# Para este script, anotaremos todos los MAGs de alta calidad generados por DAS_Tool
echo "[FASE 6] Iniciando asignación taxonómica (GTDB-Tk)..."
conda activate ont_prueba1
gtdbtk classify_wf --genome_dir $OUT_DIR/04_binning/dastool_DASTool_bins/ --out_dir $OUT_DIR/06_annotation/gtdbtk/ --extension fa --cpus $THREADS

echo "[FASE 6] Iniciando asignación taxonómica (Kraken2) y anotación funcional (Bakta)..."
conda activate anotacion_env

for mag in $OUT_DIR/04_binning/dastool_DASTool_bins/*.fa; do
    mag_name=$(basename "$mag" .fa)
    
    # Kraken2
    kraken2 --db ~/Doctorado/bases_datos/kraken2_db --threads $THREADS --use-names --report $OUT_DIR/06_annotation/kraken2_report_${mag_name}.txt $mag > $OUT_DIR/06_annotation/kraken2_${mag_name}.kraken
    
    # Bakta (Omitiendo sORFs)
    bakta --db ~/Doctorado/bases_datos/bakta_db_light/db-light --threads $THREADS --skip-sorf --output $OUT_DIR/06_annotation/bakta_${mag_name} --prefix ${mag_name} --force $mag
done

echo "[COMPLETADO] El pipeline metagenómico ha finalizado exitosamente."
