#!/usr/bin/env cwl-runner

class: CommandLineTool

cwlVersion: v1.0

baseCommand: [ "bwa", "mem" ]

requirements:
  - class: DockerRequirement
    dockerImageId: biocontainers/bwa:0.7.15
    dockerPull: biocontainers/bwa:0.7.15

inputs:
  reference_index:
    type: File
    inputBinding:
      position: 1
  fastq1_file:
    type: File
    inputBinding:
      position: 2
  fastq2_file:
    type: File
    inputBinding:
      position: 3

outputs:
  aligned_sam:
    type: stdout

stdout: aligned.sam
