#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "alignment workflow"

inputs:
  bam:
    type: File
    doc: bam file to align
  reference:
    type: File
    doc: gzipped reference fasta to align to

outputs:
  index_ref_out:
    type: File
    outputSource: index_ref/fasta_index
  bam_index_out:
    type: File
    outputSource: index_bam/bam_index

steps:
  gnu_unzip:
    run: gunzip.cwl
    in:
      reference_file: reference
    out: [ unzipped_fasta ]
  index_ref:
    run: index_fa.cwl
    in:
      reference_file: gnu_unzip/unzipped_fasta
    out: [ fasta_index ]
  sam2fastq:
    run: sam2fastq.cwl
    in:
      bam_file: bam
    out:  [ fastq1, fastq2 ]
  bwa_index:
    run: bwa_index.cwl
    in:
      reference_file: gnu_unzip/unzipped_fasta
    out: [ bwa_ref_index ]
  align_fastq:
    run: bwa_mem.cwl
    in:
      reference_index: bwa_index/bwa_ref_index
      fastq1_file: sam2fastq/fastq1
      fastq2_file: sam2fastq/fastq2
    out: [ aligned_sam ]
  sam2bam:
    run: sam2bam.cwl
    in:
      sam_file: align_fastq/aligned_sam
    out: [ aligned_bam ]
  sort_bam:
    run: sort_bam.cwl
    in:
      bam_file: sam2bam/aligned_bam
    out: [ sorted_bam ]
  index_bam:
    run: bam_index.cwl
    in:
      bam_file: sort_bam/sorted_bam
    out: [ bam_index ]
