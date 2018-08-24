#!/usr/bin/env cwl-runner

class: CommandLineTool

cwlVersion: v1.0

baseCommand: [ "gunzip" ]

arguments: [ "-c" ]

requirements:
  - class: DockerRequirement
    dockerImageId: ubuntu:xenial
    dockerPull: ubuntu:xenial

inputs:
    reference_file:
        type: File
        inputBinding:
            position: 1

outputs:
    unzipped_fasta:
        type: stdout

stdout: reference.fa
