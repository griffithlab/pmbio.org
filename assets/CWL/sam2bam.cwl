#!/usr/bin/env cwl-runner

class: CommandLineTool

cwlVersion: v1.0

baseCommand: [ "view" ]

requirements:
  - class: DockerRequirement
    dockerImageId: mgibio/samtools:1.9
    dockerPull: mgibio/samtools:1.9

arguments:
  - valueFrom: "-b"
    position: 1

inputs:
  sam_file:
    type: File
    inputBinding:
      position: 2

outputs:
  aligned_bam:
    type: stdout

stdout: aligned.bam
