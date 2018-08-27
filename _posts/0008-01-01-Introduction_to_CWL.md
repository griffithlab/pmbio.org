---
feature_text: |
  ## Precision Medicine
title: Introduction to the Common Workflow Language
categories:
    - Module 8
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-01-01
---

<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css" integrity="sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO" crossorigin="anonymous">

# Introduction to the Common Workflow Language

Often in data analysis a collection of programatic tools and resources are required to perform a specific task. This can take the form of bash scripts linking tools together however bash scripts offer no standardization and often need to be tweaked from one project to the next. The Common Workflow Language (CWL) is a specification for designing portable and scaleable workflows. It is open source and available on [github](https://github.com/common-workflow-language/common-workflow-language) under an apache 2.0 license. Using CWL instead of bash offers a number of advantages, these include:

1. Automated paralization of workflow sub tasks
2. Modular, can easily add and remove tools from a workflow
3. Portability, workflows will work with any environment with CWL installed
4. No need to monitor workflows, one sub-task failure won't cause the entire workflow to fail
5. Standardized format makes reading workflows easier

In this module we will use CWL with Docker to build an analysis pipeline to perform a simple DNA alignment.

# Installing CWL, Docker and data prerequisites

In order to begin we will need to have both CWL and docker installed. Instructions for installing both are available here:

1. [How to install docker](https://docs.docker.com/docker-for-mac/install/)
2. [How to install cwl](https://github.com/common-workflow-language/cwltool)

Further we will need reads to perform the alignment on and a reference file to align to. For this tutorial we will use downsampled reads from the HCC1395 data set. For a reference file we'll just align to chromosome 22 so things go a bit faster, we can download the reference from ensembl.

1. [HCC1395 data from a single lane](https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395_1tenth_percent/gerald_C1TD1ACXX_7_ATCACG.bam)
2. [Chromosome 22 Fasta](ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz)

# CWL Workflow Pieces

A typical CWL workflow consists of three main pieces, the first piece is simply a yaml file specifying the inputs to the workflow. In this tutorial the inputs are simply the bam file containing reads to align and a fasta file to align the reads to. The second piece is a cwl file containing the workflow, in other words how things are run, what the outputs should be, etc. The last piece are cwl files specifying how the tools will be run. Don't worry if this all doesn't make sense, things should clear up as we go along. For our example as mentioned we will be constructing a workflow to perform DNA alignment. Go ahead and download the yml and cwl files and put them all in the same directory. You can do so by clicking on the links below:

1. data inputs
    1. [gerald_C1TD1ACXX_7_ATCACG.bam](https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395_1tenth_percent/gerald_C1TD1ACXX_7_ATCACG.bam)
    2. [Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz](ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz)

2. yaml file specifying the inputs
    1. [inputs.yml](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/inputs.yml)

3. Workflow file
    1. [workflow.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/workflow.cwl)

4. Individual files to run tools
    1. [gunzip.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/gunzip.cwl)
    2. [index_fa.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/index_fa.cwl)
    3. [sam2fastq.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/sam2fastq.cwl)
    4. [bwa_index.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/bwa_index.cwl)
    5. [bwa_mem.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/bwa_mem.cwl)
    6. [sam2bam.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/sam2bam.cwl)
    7. [sort_bam.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/sort_bam.cwl)
    8. [bam_index.cwl](https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/CWL/bam_index.cwl)

# The inputs.yml file

Okay lets start by going over what the input.yml file is. Simply put, as it sounds it is simply specifying the inputs given to the workflow. In our workflow we only have two inputs, a bam file and a reference file. The inputs.yml is specifying what the input is (i.e. files), where the inputs exist (i.e. file paths) and the identifier the cwl workflow will use to refer to the inputs.

<div class="highlight"><pre class="highlight"><code><span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="Identifier for input bam file">bam</span><span class="pi">:</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="Field indicating the object type to expect">class</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates that input is a file">File</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating where the file is located">path</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="path to the input bam file">/Users/zskidmor/Desktop/lab_meeting/gerald_C1TD1ACXX_7_ATCACG.bam</span>
<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="Identifier for input reference file">reference</span><span class="pi">:</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="Field indicating the object type to expect">class</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the input is a file">File</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating where the file is located">path</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="path to the fasta file">/Users/zskidmor/Desktop/lab_meeting/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz</span>
</code></pre></div>

# The workflow.cwl files

The workflow file specifies how things should be run, the inputs, outputs, and steps corresponding to the specific workflow.

<div class="language-yaml highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="c1" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="shebang line indicating the cwl runner to use" >#!/usr/bin/env cwl-runner</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the cwl version, this is important as small differences in available fields, etc. could be different between different versions" >cwlVersion</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates cwl version is 1.0">v1.0</span>
<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the workflow type, one of Workflow or Command Line Tool">class</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="Indicates that this is a workflow with steps, etc.">Workflow</span>
<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="label for the workflow being run">label</span><span class="pi">:</span> <span class="s2">"</span><span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates that this is a workflow for alignment">alignment</span><span class="nv"> </span><span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates this is a workflow for alignment">workflow"</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the inputs to the workflow">inputs</span><span class="pi">:</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="refers to the bam identifier from the input.yml file above">bam</span><span class="pi">:</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the type of object to expect">type</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the type of object to expect is a file">File</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field referencing documentation for the bam input">doc</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="documents the bam input">bam file to align</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="refers to the reference identifier from the input.yml file above">reference</span><span class="pi">:</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the type of object to expect">type</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the type of object to expect is a file">File</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field referencing documentation for the reference input">doc</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="documents the reference input">gzipped reference fasta to align to</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the outputs to the workflow">outputs</span><span class="pi">:</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="identifier for the output for the reference file">index_ref_out</span><span class="pi">:</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating what the output type is">type</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates that the output is a file">File</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field for where the output is comming from">outputSource</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the output is coming from the index_ref step and the output is the fasta_index from that step">index_ref/fasta_index</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="identifier for the output for the bam file">bam_index_out</span><span class="pi">:</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating what the output type is">type</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the output is a file">File</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field for where the output is comming from">outputSource</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the output is coming from the index_bam step and the output is the bam_index from that step">index_bam/bam_index</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates steps to the workflow">steps</span><span class="pi">:</span>
  <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="name for gunzip step">gnu_unzip</span><span class="pi">:</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field specifying what should be run">run</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the cwl file to be run in this step">gunzip.cwl</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field specifying input for this step">in</span><span class="pi">:</span>
      <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the identifier for the input for this step">reference_file</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the input for this step is the reference identifier from the inputs section">reference</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field specifying the output identifier from this step">out</span><span class="pi">:</span> <span class="pi">[</span> <span class="nv" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the identifier for the output which is referred to in gunzip.cwl">unzipped_fasta</span> <span class="pi">]</span>
  <span class="na">index_ref</span><span class="pi">:</span>
    <span class="na">run</span><span class="pi">:</span> <span class="s">index_fa.cwl</span>
    <span class="na">in</span><span class="pi">:</span>
      <span class="na">reference_file</span><span class="pi">:</span> <span class="s">gnu_unzip/unzipped_fasta</span>
    <span class="na">out</span><span class="pi">:</span> <span class="pi">[</span> <span class="nv">fasta_index</span> <span class="pi">]</span>
  <span class="na">sam2fastq</span><span class="pi">:</span>
    <span class="na">run</span><span class="pi">:</span> <span class="s">sam2fastq.cwl</span>
    <span class="na">in</span><span class="pi">:</span>
      <span class="na">bam_file</span><span class="pi">:</span> <span class="s">bam</span>
    <span class="na">out</span><span class="pi">:</span>  <span class="pi">[</span> <span class="nv">fastq1</span><span class="pi">,</span> <span class="nv">fastq2</span> <span class="pi">]</span>
  <span class="na">bwa_index</span><span class="pi">:</span>
    <span class="na">run</span><span class="pi">:</span> <span class="s">bwa_index.cwl</span>
    <span class="na">in</span><span class="pi">:</span>
      <span class="na">reference_file</span><span class="pi">:</span> <span class="s">gnu_unzip/unzipped_fasta</span>
    <span class="na">out</span><span class="pi">:</span> <span class="pi">[</span> <span class="nv">bwa_ref_index</span> <span class="pi">]</span>
  <span class="na">align_fastq</span><span class="pi">:</span>
    <span class="na">run</span><span class="pi">:</span> <span class="s">bwa_mem.cwl</span>
    <span class="na">in</span><span class="pi">:</span>
      <span class="na">reference_index</span><span class="pi">:</span> <span class="s">bwa_index/bwa_ref_index</span>
      <span class="na">fastq1_file</span><span class="pi">:</span> <span class="s">sam2fastq/fastq1</span>
      <span class="na">fastq2_file</span><span class="pi">:</span> <span class="s">sam2fastq/fastq2</span>
    <span class="na">out</span><span class="pi">:</span> <span class="pi">[</span> <span class="nv">aligned_sam</span> <span class="pi">]</span>
  <span class="na">sam2bam</span><span class="pi">:</span>
    <span class="na">run</span><span class="pi">:</span> <span class="s">sam2bam.cwl</span>
    <span class="na">in</span><span class="pi">:</span>
      <span class="na">sam_file</span><span class="pi">:</span> <span class="s">align_fastq/aligned_sam</span>
    <span class="na">out</span><span class="pi">:</span> <span class="pi">[</span> <span class="nv">aligned_bam</span> <span class="pi">]</span>
  <span class="na">sort_bam</span><span class="pi">:</span>
    <span class="na">run</span><span class="pi">:</span> <span class="s">sort_bam.cwl</span>
    <span class="na">in</span><span class="pi">:</span>
      <span class="na">bam_file</span><span class="pi">:</span> <span class="s">sam2bam/aligned_bam</span>
    <span class="na">out</span><span class="pi">:</span> <span class="pi">[</span> <span class="nv">sorted_bam</span> <span class="pi">]</span>
  <span class="na">index_bam</span><span class="pi">:</span>
    <span class="na">run</span><span class="pi">:</span> <span class="s">bam_index.cwl</span>
    <span class="na">in</span><span class="pi">:</span>
      <span class="na">bam_file</span><span class="pi">:</span> <span class="s">sort_bam/sorted_bam</span>
    <span class="na">out</span><span class="pi">:</span> <span class="pi">[</span> <span class="nv">bam_index</span> <span class="pi">]</span>
</code></pre></div></div>

# Command.yml files

The command.yml files specify how to run a given command for a step in the workflow. In the example below we go over how the file is structured for the gnu_unzip step specified in the workflow.yaml.

<div class="language-yaml highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="c1" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="shebang line indicating the runner">#!/usr/bin/env cwl-runner</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the type of cwl file">class</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="specifies that this cwl file is for running a command not a workflow">CommandLineTool</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating what cwl version file corresponds to">cwlVersion</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates that the cwl version in the file corresponds to v1.0">v1.0</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the base command to run">baseCommand</span><span class="pi">:</span> <span class="pi">[</span> <span class="s2">"</span><span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the base command in gunzip">gunzip"</span> <span class="pi">]</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="specifies any arguments to be added to the base command">arguments</span><span class="pi">:</span> <span class="pi">[</span> <span class="s2">"</span><span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="specifies gunzip should be run with the -c option which indicates it should print to stdout">-c"</span> <span class="pi">]</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating any additional requirements when running the gunzip command">requirements</span><span class="pi">:</span>
  <span class="pi">-</span> <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the class of the requirement">class</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates that the command needs to run with docker">DockerRequirement</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the name of the docker image on your computer">dockerImageId</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the docker image name is ubuntu:xenial">ubuntu:xenial</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the name of the docker image to pull if no image id is found locally">dockerPull</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the docker image name to pull is ubuntu:xenial">ubuntu:xenial</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field specifying the inputs for this command">inputs</span><span class="pi">:</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the identifier associated with the input given in workflow.cwl">reference_file</span><span class="pi">:</span>
        <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field specifying the class of input to expect">type</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the input is expected to be a file">File</span>
        <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field specifying various options for the input">inputBinding</span><span class="pi">:</span>
            <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the position of this input relative to the base command">position</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates this input is at position 1 i.e. right after the base command and arguments">1</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field specifying the output from this command">outputs</span><span class="pi">:</span>
    <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="identifier for the output of this command">unzipped_fasta</span><span class="pi">:</span>
        <span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field indicating the type of output to expect">type</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates the type of output to pbe printed to stdout">stdout</span>

<span class="na" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="field for capturing output in stdout">stdout</span><span class="pi">:</span> <span class="s" tabindex="0" data-toggle="popover" data-trigger="focus" data-content="indicates that stdout should go to a file named reference.fa">reference.fa</span>
</code></pre></div></div>

<script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.3/umd/popper.min.js" integrity="sha384-ZMP7rVo3mIykV+2+9J3UJ46jBk0WLaUAdn689aCwoqbBJiSnjAK/l8WvCWPIPm49" crossorigin="anonymous"></script>
<script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/js/bootstrap.min.js" integrity="sha384-ChfqqxuZUCnJSK3+MXmPNIyE6ZbWh2IMqE241rYiqJxyMiZ6OW/JmZQ5stwEULTy" crossorigin="anonymous"></script>

# Putting it all together

Now that we've gone over the basics let's go ahead and run this workflow. On a typical computer the workflow should run in aprox. 7-10 minutes depending on if docker images need to be pulled down from the web. 

```bash
cwltool --outdir ~/Desktop/cwl_test workflow.cwl inputs.yml
```

<script>
$(function () {
  $('[data-toggle="popover"]').popover()
})
</script>

<script>
$('.popover-dismiss').popover({
  trigger: 'focus'
})
</script>
