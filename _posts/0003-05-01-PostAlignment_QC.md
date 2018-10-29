---
feature_text: |
  ## Precision Medicine
title: Post-Alignment QC
categories:
    - Module-03-Align
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-05-01
---

Calculate QC metric for bam files using samtools and picard

### Run Samtools flagstat

```bash
samtools flagstat /workspace/align/Exome_Norm_merged_sorted_mrkdup_bqsr.bam > /workspace/align/Exome_Norm_flagstat.txt
samtools flagstat /workspace/align/${TUMOR_DATA_SM}_exome.bam > /workspace/align${TUMOR_DATA_SM}_exome_flagstat.txt
samtools flagstat /workspace/align/${NORMAL_DATA_SM}_wgs.bam > /workspace/align/${NORMAL_DATA_SM}_wgs_flagstat.txt
samtools flagstat /workspace/align/${TUMOR_DATA_SM}_wgs.bam > /workspace/align/${TUMOR_DATA_SM}_wgs_flagstat.txt
```
### Run various picard CollectMetrics tools

```bash
java -Xmx24g -jar $PICARD CollectInsertSizeMetrics I=/workspace/align/${NORMAL_DATA_SM}_exome.bam O=/workspace/align/${NORMAL_DATA_SM}_exome_insert_size_metrics.txt H=/workspace/align/${NORMAL_DATA_SM}_exome_insert_size_metrics.pdf
$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectInsertSizeMetrics I=/workspace/align/${TUMOR_DATA_SM}_exome.bam O=/workspace/align/${TUMOR_DATA_SM}_exome_insert_size_metrics.txt H=/workspace/align/${TUMOR_DATA_SM}_exome_insert_size_metrics.pdf
$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectInsertSizeMetrics I=/workspace/align/${NORMAL_DATA_SM}_wgs.bam O=/workspace/align/${NORMAL_DATA_SM}_wgs_insert_size_metrics.txt H=/workspace/align/${NORMAL_DATA_SM}_wgs_insert_size_metrics.pdf
$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectInsertSizeMetrics I=/workspace/align/${TUMOR_DATA_SM}_wgs.bam O=/workspace/align/${TUMOR_DATA_SM}_wgs_insert_size_metrics.txt H=/workspace/align/${TUMOR_DATA_SM}_wgs_insert_size_metrics.pdf

# Picard CollectAlignmentSummaryMetrics
$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectAlignmentSummaryMetrics I=/workspace/align/${NORMAL_DATA_SM}_exome.bam O=/workspace/align/${NORMAL_DATA_SM}_exome_alignment_metrics.txt R=$SOMATIC_REFSEQ_FASTA
$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectAlignmentSummaryMetrics I=/workspace/align/${TUMOR_DATA_SM}_exome.bam O=/workspace/align/${TUMOR_DATA_SM}_exome_alignment_metrics.txt R=$SOMATIC_REFSEQ_FASTA
$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectAlignmentSummaryMetrics I=/workspace/align/${NORMAL_DATA_SM}_wgs.bam O=/workspace/align/${NORMAL_DATA_SM}_wgs_alignment_metrics.txt R=$SOMATIC_REFSEQ_FASTA
$JAVA_EIGHT -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectAlignmentSummaryMetrics I=/workspace/align/${TUMOR_DATA_SM}_wgs.bam O=/workspace/align/${TUMOR_DATA_SM}_wgs_alignment_metrics.txt R=$SOMATIC_REFSEQ_FASTA

#Picard CollectHsMetrics
#Exome Only
$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectHsMetrics I=/workspace/align/${NORMAL_DATA_SM}_exome.bam O=/workspace/align/${NORMAL_DATA_SM}_exome_hs_metrics.txt R=$SOMATIC_REFSEQ_FASTA BI=$SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-probes.interval_list TI=$SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-targets.interval_list
$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectHsMetrics I=/workspace/align/${TUMOR_DATA_SM}_exome.bam O=/workspace/align/${TUMOR_DATA_SM}_exome_hs_metrics.txt R=$SOMATIC_REFSEQ_FASTA BI=$SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-probes.interval_list TI=$SOMATIC_REFSEQ_DIR/xgen-exome-research-panel-targets.interval_list

#Picard CollectGcBiasMetrics
#WGS only
$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectGcBiasMetrics I=/workspace/align/${NORMAL_DATA_SM}_wgs.bam O=/workspace/align/${NORMAL_DATA_SM}_wgs_gc_bias_metrics.txt R=$SOMATIC_REFSEQ_FASTA CHART=/workspace/align/${NORMAL_DATA_SM}_wgs_gc_bias_metrics.pdf S=/workspace/align/${NORMAL_DATA_SM}_wgs_gc_bias_summary.txt
$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectGcBiasMetrics I=/workspace/align/${TUMOR_DATA_SM}_wgs.bam O=/workspace/align/${TUMOR_DATA_SM}_wgs_gc_bias_metrics.txt R=$SOMATIC_REFSEQ_FASTA CHART=/workspace/align/${TUMOR_DATA_SM}_wgs_gc_bias_metrics.pdf S=/workspace/align/${TUMOR_DATA_SM}_wgs_gc_bias_summary.txt

#Picard CollectWgsMetrics
$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectWgsMetrics I=/workspace/align/${NORMAL_DATA_SM}_wgs.bam O=/workspace/align/${NORMAL_DATA_SM}_wgs_metrics.txt R=$SOMATIC_REFSEQ_FASTA INTERVALS=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_autosomal.interval_list
$JAVA_EIGHT -Xmx${JVM_RAM_GB}g -jar $SOMATIC_HOME/software/picard-tools-2.4.1/picard.jar CollectWgsMetrics I=/workspace/align/${TUMOR_DATA_SM}_wgs.bam O=/workspace/align/${TUMOR_DATA_SM}_wgs_metrics.txt R=$SOMATIC_REFSEQ_FASTA INTERVALS=$SOMATIC_REFSEQ_DIR/${GENOME_BASENAME}_autosomal.interval_list
```
### Run FastQC


### Run MultiQC to produce a final report


See docs here: https://github.com/genome/cancer-genomics-workflow/wiki/Alignment
