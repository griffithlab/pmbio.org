#!/usr/bin/env cwl-runner

class: CommandLineTool

cwlVersion: v1.0

baseCommand: [ "bwa", "index" ]

inputs:
  reference_file:
    type: File
    inputBinding:
      position: 1

requirements:
  - class: DockerRequirement
    dockerPull: biocontainers/bwa:0.7.15
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.reference_file)

outputs:
  bwa_ref_index:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
    outputBinding:
      glob: "*.fa"
