#!/usr/bin/env cwl-runner

class: CommandLineTool

cwlVersion: v1.0

baseCommand: [ "faidx" ]

requirements:
  - class: DockerRequirement
    dockerPull: mgibio/samtools:1.9
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.reference_file)

inputs:
  reference_file:
    type: File
    inputBinding:
      position: 1

outputs:
  fasta_index:
    type: File
    secondaryFiles: [.fai]
    outputBinding:
      glob: "*.fa"
