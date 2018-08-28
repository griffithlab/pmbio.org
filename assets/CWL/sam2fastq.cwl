#!/usr/bin/env cwl-runner

class: CommandLineTool

cwlVersion: v1.0

baseCommand: [ "SamToFastq" ]

arguments:
- valueFrom: reads.1.fastq
  position: 2
  prefix: FASTQ=
- valueFrom: reads.2.fastq
  position: 3
  prefix: SECOND_END_FASTQ=

requirements:
  - class: DockerRequirement
    dockerPull: broadinstitute/picard:2.18.2

inputs:
    bam_file:
        type: File
        inputBinding:
            position: 1
            prefix: I=

outputs:
  fastq1:
    type: File
    outputBinding:
      glob: reads.1.fastq
  fastq2:
    type: File
    outputBinding:
      glob: reads.2.fastq
