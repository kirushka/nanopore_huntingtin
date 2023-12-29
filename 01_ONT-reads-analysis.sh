#!/bin/bash

################################################################################
# Combine all reads into one FASTQ file
cat fastq_pass/FAR91602_pass_barcode04_9625fc3e_0.fastq fastq_pass/FAR91602_pass_barcode05_9625fc3e_0.fastq fastq_pass/FAR91602_pass_barcode06_9625fc3e_0.fastq fastq_pass/FAR91602_pass_barcode07_9625fc3e_0.fastq fastq_pass/FAR91602_pass_barcode08_9625fc3e_0.fastq fastq_pass/FAR91602_pass_unclassified_9625fc3e_0.fastq > fastq_pass/fastq_pass.fastq


################################################################################
# QC with NanoPlot
conda activate nanoplot 
NanoPlot -t 6 --fastq ../../fastq_pass.fastq --raw --plots dot --tsv_stats -p ALL -o fastq_pass/QC/nanoplot_raw_barcoded


################################################################################
# Group barcoded samples into Q15 or Q138 group: 
# Q15 - barcodes 7,8 
cat fastq_pass/FAR91602_pass_barcode07_9625fc3e_0.fastq fastq_pass/FAR91602_pass_barcode08_9625fc3e_0.fastq > fastq_pass/combined_barcodes/Q15_combined.fastq

# Q138 - barcodes 4,5,6
cat fastq_pass/FAR91602_pass_barcode04_9625fc3e_0.fastq fastq_pass/FAR91602_pass_barcode05_9625fc3e_0.fastq fastq_pass/FAR91602_pass_barcode06_9625fc3e_0.fastq > fastq_pass/combined_barcodes/Q138_combined.fastq

# Keep unclassified reads as an individual group
cp fastq_pass/FAR91602_pass_unclassified_9625fc3e_0.fastq fastq_pass/combined_barcodes/UNC.fastq


################################################################################
# Align reads to linearized reference plasmid sequence
conda activate minimap2

# Q15 reads
minimap2 -ax map-ont reference/pSBtet-Neo-HttQ15_linear.fa fastq_pass/combined_barcodes/Q15_combined.fastq | samtools view -bS - | samtools sort -o align_minimap2/combined_barcodes/Q15_plasmid_linear.psort.bam
samtools index align_minimap2/combined_barcodes/Q15_plasmid_linear.psort.bam

# Q138 reads
minimap2 -ax map-ont reference/pSBtet-Neo-HttQ15_linear.fa fastq_pass/combined_barcodes/Q138_combined.fastq | samtools view -bS - | samtools sort -o align_minimap2/combined_barcodes/Q138_plasmid_linear.psort.bam
samtools index align_minimap2/combined_barcodes/Q38_plasmid_linear.psort.bam

# Unclassified reads
minimap2 -ax map-ont reference/pSBtet-Neo-HttQ15_linear.fa fastq_pass/combined_barcodes/UNC.fastq | samtools view -bS - | samtools sort -o align_minimap2/combined_barcodes/UNC_plasmid_linear.psort.bam
samtools index align_minimap2/combined_barcodes/UNC_plasmid_linear.psort.bam


################################################################################
# Keep reads that align to N-terminal end of Htt gene (position 1-147) and convert BAM to FASTQ
# Q15 reads
bedtools intersect -wa -abam align_minimap2/combined_barcodes/Q15_plasmid_linear.psort.bam -b reference/N_terminal_htt_in_plasmid.bed -F 1.0 > align_minimap2/combined_barcodes/Q15_plasmid_linear.psort.flt.bam
samtools index align_minimap2/combined_barcodes/Q15_plasmid_linear.psort.flt.bam
samtools fastq -F 0 align_minimap2/combined_barcodes/Q15_plasmid_linear.psort.flt.bam > align_minimap2/combined_barcodes/fastq/Q15_plasmid_linear.flt.fastq

# Q138 reads
bedtools intersect -wa -abam align_minimap2/combined_barcodes/Q138_plasmid_linear.psort.bam -b reference/N_terminal_htt_in_plasmid.bed -F 1.0 > align_minimap2/combined_barcodes/Q138_plasmid_linear.psort.flt.bam
samtools index align_minimap2/combined_barcodes/Q138_plasmid_linear.psort.flt.bam
samtools fastq -F 0 align_minimap2/combined_barcodes/Q138_plasmid_linear.psort.flt.bam > align_minimap2/combined_barcodes/fastq/Q138_plasmid_linear.flt.fastq

# Unclassified reads
bedtools intersect -wa -abam align_minimap2/combined_barcodes/UNC_plasmid_linear.psort.bam -b reference/N_terminal_htt_in_plasmid.bed -F 1.0 > align_minimap2/combined_barcodes/UNC_plasmid_linear.psort.flt.bam
samtools index align_minimap2/combined_barcodes/UNC_plasmid_linear.psort.flt.bam
samtools fastq -F 0 align_minimap2/combined_barcodes/UNC_plasmid_linear.psort.flt.bam > align_minimap2/combined_barcodes/fastq/UNC_plasmid_linear.flt.fastq

# Combine all filtered reads into one FASTQ file
cat align_minimap2/combined_barcodes/fastq/Q15_plasmid_linear.flt.fastq align_minimap2/combined_barcodes/fastq/Q138_plasmid_linear.flt.fastq align_minimap2/combined_barcodes/fastq/UNC_plasmid_linear.flt.fastq > align_minimap2/combined_barcodes/fastq/ALL_plasmid_linear.flt.fastq


################################################################################
# CAG repeat length ensimation
conda activate last_tandem_genotypes

# Index reference
lastdb -uNEAR reference/htt_linear reference/pSBtet-Neo-HttQ15_linear.fa

# Calculate substitution and gap rates
last-train -Q0  reference/htt_linear align_minimap2/combined_barcodes/fastq/ALL_plasmid_linear.flt.fastq > align_last/plasmid_flt/plasmid_linear.par

# Align reads to reference
# Q15 reads
lastal --split -p align_last/plasmid_flt/plasmid_linear.par reference/htt_linear align_minimap2/combined_barcodes/fastq/Q15_plasmid_linear.flt.fastq > align_last/plasmid_flt/Q15_plasmid_linear.maf
# Q138 reads
lastal --split -p align_last/plasmid_flt/plasmid_linear.par reference/htt_linear align_minimap2/combined_barcodes/fastq/Q138_plasmid_linear.flt.fastq > align_last/plasmid_flt/Q138_plasmid_linear.maf
# Unclassified reads
lastal --split -p align_last/plasmid_flt/plasmid_linear.par reference/htt_linear align_minimap2/combined_barcodes/fastq/UNC_plasmid_linear.flt.fastq > align_last/plasmid_flt/UNC_plasmid_linear.maf

# Run tandem-genotypes (it requires bed-like file with repeat annotation)
tandem-genotypes htt_repeat_linear.txt align_last/plasmid_flt/Q15_plasmid_linear.maf > tandem_genotypes/plasmid_flt/Q15_plasmid_linear.txt
tandem-genotypes htt_repeat_linear.txt align_last/plasmid_flt/Q138_plasmid_linear.maf > tandem_genotypes/plasmid_flt/Q138_plasmid_linear.txt
tandem-genotypes htt_repeat_linear.txt align_last/plasmid_flt/UNC_plasmid_linear.maf > tandem_genotypes/plasmid_flt/UNC_plasmid_linear.txt
