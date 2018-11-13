---
feature_text: |
  ## Precision Medicine Bioinformatics
  Introduction to bioinformatics for DNA and RNA sequence analysis
title: Intro to IGV
categories:
    - Module-03-Align
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-01
---

It is often necessary to examine sequencing data aligned to specific regions of the genome in order to obtain a clearer picture of genomic events. One of the most popular tools for this is the [Integrative Genomics Viewer](http://software.broadinstitute.org/software/igv/). After this lab you should be able to perform the following tasks:
1. Visualize a variety of genomic data
2. Quickly navigate around the genome
3. Visualize read alignments
4. Validate SNP/SNV calls and structural re-arrangements by eye

### Installing IGV

Java is necessary to run IGV, you can download the java runtime environment (JRE) for your operating system [here](http://www.oracle.com/technetwork/java/javase/downloads/index.html). To determine if this step is necessary type `java -version` at a command prompt, if the program is not >= 1.7 you'll need to upgrade it. IGV can be downloaded [here](https://software.broadinstitute.org/software/igv/download). This tutorial will make use of IGV version 2.3, we strongly recommend that you upgrade IGV if you have an older version installed.

### Data Set for IGV
We will be using publicly available Illumina sequence data from the HCC1143 cell line. The HCC1143 cell line was generated from a 52 year old caucasian woman with breast cancer. Additional information on this cell line can be found [here](http://www.atcc.org/products/all/CRL-2321.aspx): (tumor, TNM stage IIA, grade 3, primary ductal carcinoma) and [HCC1143/BL](http://www.atcc.org/products/all/CRL-2362.aspx) (matched normal EBV transformed lymphoblast cell line). Reads within these cell lines have been filtered to Chromosome 21: 19,000,000-20,000,000 in order to reduce file sizes.

* [HCC1143.normal.21.19M-20M.bam](http://genomedata.org/gen-viz-workshop/IGV/HCC1143.normal.21.19M-20M.bam)
* [HCC1143.normal.21.19M-20M.bam.bai](http://genomedata.org/gen-viz-workshop/IGV/HCC1143.normal.21.19M-20M.bam.bai)

### Visualization Part 1: Getting familiar with IGV

We will be visualizing read alignments using [IGV](https://software.broadinstitute.org/software/igv/), a popular visualization tool for HTS data.

First, lets familiarize ourselves with it.

#### Load a Genome and some Data Tracks

By default, IGV loads Human (hg19). If you work with another version of the human genome, or another organism altogether, you can change the genome by clicking the drop down menu in the upper-left. For this lab, we will be using Human (hg19).

We will also load additional tracks from the IGV **Server** using (`File` -> `Load from Server...`):

* Ensembl genes (or your favourite source of gene annotations)
* GC Percentage
* dbSNP 1.3.1 or 1.3.7

**Load hg19 genome and additional data tracks**
{% include figure.html image="/assets/IGV/load.data.tracks.png" width="450" %}

#### Navigation

You should see a listing of chromosomes for this reference genome. Choose ***1***, for chromosome 1.

**Chromosome chooser**
{% include figure.html image="/assets/IGV/chromosomes.png" width="750" %}


Navigate to **chr1:10,000-11,000** by entering this into the location field (in the top-left corner of the interface) and clicking `Go`. This shows a window of chromosome 1 that is 1,000 base pairs wide and beginning at position 10,000.

**Navigition using Location text field. Sequence displayed as thin coloured rectangles.**
{% include figure.html image="/assets/IGV/1.png" width="750" %}


IGV displays the sequence of letters in a genome as a sequence of colours (e.g. A = green, C = blue, etc.). This makes repetitive sequences, like the ones found at the start of this region, easy to identify. Zoom in a bit more using the `+` button (top right) to see the individual bases of the reference genome sequence.

You can navigate to a gene of interest by typing it into the same box that the genomic coordinates are in and pressing Enter/Return. Try it for your favourite gene, or *BRCA1* if you can not decide.

**Gene model**
{% include figure.html image="/assets/IGV/gene_model.png" width="450" %}

Genes are represented as lines and boxes. Lines represent intronic regions, and boxes represent exonic regions. The arrows indicate the direction/strand of transcription for the gene. When an exon box become narrower in height, this indicates a UTR.

When loaded, tracks are stacked on top of each other. You can identify which track is which by consulting the label to the left of each track.

#### Region Lists

Sometimes, it is useful to save where you are, or to load regions of interest. For this purpose, there is a **Region Navigator** in IGV. To access it, click `Regions` > `Region Navigator`. While you browse around the genome, you can save some bookmarks by pressing the `Add` button at any time.

**Bookmarks in IGV**
{% include figure.html image="/assets/IGV/bookmarks.png" width="450" %}

#### Loading Read Alignments
We will be using the breast cancer cell line HCC1143 to visualize alignments. For speed, only a small portion of chr21 will be loaded (19M:20M).

**HCC1143 Alignments to hg19:**
* [HCC1143.normal.21.19M-20M.bam](http://genomedata.org/gen-viz-workshop/IGV/HCC1143.normal.21.19M-20M.bam)
* [HCC1143.normal.21.19M-20M.bam.bai](http://genomedata.org/gen-viz-workshop/IGV/hcc1143/HCC1143.normal.21.19M-20M.bam.bai)

Copy the files to your local drive, and in IGV choose `File` > `Load from File...`, select the bam file, and click `OK`. Note that the bam and index files must be in the same directory for IGV to load these properly. Alternatively, you can copy the link location and load `File` > `Load from URL...`.

**Load BAM track from File**
{% include figure.html image="/assets/IGV/load_bam.png" width="750" %}

#### Visualizing read alignments

Navigate to a narrow window on chromosome 21: `chr21:19,480,041-19,480,386`.

To start our exploration, right click on the track-name, and select the following options:
* Sort alignments by `start location`
* Group alignments by `pair orientation`

Experiment with the various settings by right clicking the read alignment track and toggling the options. Think about which would be best for specific tasks (e.g. quality control, SNP calling, CNV finding).

**Changing how read alignments are sorted, grouped, and colored**<br>
{% include figure.html image="/assets/IGV/sort_and_group.png" width="450" %}

You will see reads represented by grey or white bars stacked on top of each other, where they were aligned to the reference genome. The reads are pointed to indicate their orientation (i.e. the strand on which they are mapped). Mouse over any read and notice that a lot of information is available. To toggle read display from `hover` to `click`, select the yellow box and change the setting.

**Changing how read information is shown (i.e. on hover, click, never)**
{% include figure.html image="/assets/IGV/show_details_on_click.png" width="750" %}

Once you select a read, you will learn what many of these metrics mean, and how to use them to assess the quality of your datasets.  At each base that the read sequence **mismatches** the reference, the colour of the base represents the letter that exists in the read (using the same colour legend used for displaying the reference).

**Viewing read information for a single aligned read**
{% include figure.html image="/assets/IGV/click_read.png" width="750" %}

### Visualization Part 2: Inspecting SNPs, SNVs, and SVs

In this section we will be looking in detail at 8 positions in the genome, and determining whether they represent real events or artifacts.

#### Two neighbouring SNPs

* Navigate to region `chr21:19,479,237-19,479,814`
* Note two heterozygous variants, one corresponds to a known dbSNP (`G/T` on the right) the other does not (`C/T` on the left)
* Zoom in and center on the `C/T` SNV on the left (`chr21:19,479,321` is the SNV position)
* Sort alignments by `base`
* Color alignments by `read strand`

**Example1. Good quality SNVs/SNPs**
{% include figure.html image="/assets/IGV/example1_color.png" width="750" %}

**Notes:**
* High base qualities in all reads except one (where the alt allele is the last base of the read)
* Good mapping quality of reads, no strand bias, allele frequency consistent with heterozygous mutation

{% include question.html question="What does Shade base by quality do?" answer='This will change the opacity of the base in IGV based on how confident the sequencer was in calling that base using the phred score. This is beneficial in determining if a called variant is real or artifactual.'%}
{% include question.html question="How does Color by read strand help?" answer='Coloring by read strand will indicate if the DNA fragment sequenced was on the positive or negative strand. A variant occurring on only one strand could indicate an artifact.'%}

#### Homopolymer region with indel

Navigate to position `chr21:19,518,412-19,518,497`

**Example 2a**
* Group alignments by `read strand`
* Center on the `A` within the homopolymer run (`chr21:19,518,470`), and `Sort alignments by` -> `base`

{% include figure.html image="/assets/IGV/example2a.png" width="750" %}


**Example 2b**
* Center on the one base deletion (`chr21:19,518,452`), and `Sort alignments by` -> `base`

{% include figure.html image="/assets/IGV/example2b.png" width="750" %}

**Notes:**
* The alt allele is either a deletion or insertion of one or two `T`s
* The remaining bases are mismatched, because the alignment is now out of sync
* The dpSNP entry at this location (rs74604068) is an A->T, and in all likelihood an artifact
* i.e. the common variants from dbSNP include some cases that are actually common misalignments caused by repeats

## Coverage by GC

Navigate to position `chr21:19,611,925-19,631,555`. Note that the range contains areas where coverage drops to zero in a few places.

**Example 3**
* Use `Collapsed` view
* Use `Color alignments by` -> `insert size and pair orientation`
* Load GC track
* See concordance of coverage with GC content

{% include figure.html image="/assets/IGV/example3.png" width="750" %}

{% include question.html question="Why are there blue and red reads throughout the alignments?" answer='The reads are colored by insert size, in paired data a blue read indicates the insert size is smaller than expected indicating a deletion. Conversely a red read indicates the insert size is larger than expected indicating an insertion.'%}

## Heterozygous SNPs on different alleles

Navigate to region `chr21:19,666,833-19,667,007`

**Example 4**
* Sort by base (at position `chr21:19,666,901`)

{% include figure.html image="/assets/IGV/example4.png" width="750" %}

**Note:**
* There is no linkage between alleles for these two SNPs because reads covering both only contain one or the other

#### Low mapping quality

Navigate to region `chr21:19,800,320-19,818,162`
* Load repeat track (`File` -> `Load from server...`)

**Load repeats**<br>
{% include figure.html image="/assets/IGV/load_repeats.png" width="450" %}

**Example 5**
{% include figure.html image="/assets/IGV/example5.png" width="750" %}

**Notes:**
* Mapping quality plunges in all reads (white instead of grey).  Once we load repeat elements, we see that there are two LINE elements that cause this.

#### Homozygous deletion

Navigate to region `chr21:19,324,469-19,331,468`

**Example 6**
* Turn on `View as Pairs` and `Expanded` view
* Use `Color alignments by` -> `insert size and pair orientation`
* Sort reads by insert size
* Click on a red read pair to pull up information on alignments

{% include figure.html image="/assets/IGV/example6.png" width="750" %}

**Notes:**
* Typical insert size of read pair in the vicinity: 350bp
* Insert size of red read pairs: 2,875bp
* This corresponds to a homozygous deletion of 2.5kb

#### Mis-alignment

Navigate to region `chr21:19,102,154-19,103,108`

**Example 7**
{% include figure.html image="/assets/IGV/example7.png" width="750" %}

**Notes:**
* This is a position where AluY element causes mis-alignment.
* Misaligned reads have mismatches to the reference and well-aligned reads have partners on other chromosomes where additional ALuY elements are encoded.
* Zoom out until you can clearly see the contrast between the difficult alignment region (corresponding to an AluY) and regions with clean alignments on either side

#### Translocation

Navigate to region `chr21:19,089,694-19,095,362`

**Example 8**
* Expanded view
* `Group alignments by` -> `pair orientation`
* `Color alignments by` -> `pair orientation`

{% include figure.html image="/assets/IGV/example8.png" width="750" %}

**Notes:**
* Many reads with mismatches to reference
* Read pairs in RL pattern (instead of LR pattern)
* Region is flanked by reads with poor mapping quality (white instead of grey)
* Presence of reads with pairs on other chromosomes (coloured reads at the bottom when scrolling down)

### Visualization Part 3: Automating Tasks in IGV

We can use the Tools menu to invoke running a batch script. Batch scripts are described on the IGV website:
* Batch file requirements: [http://software.broadinstitute.org/software/igv/batch](http://software.broadinstitute.org/software/igv/batch)
* Commands recognized in a batch script: [http://software.broadinstitute.org/software/igv/PortCommands](http://software.broadinstitute.org/software/igv/PortCommands)
* We also need to provide sample attribute file as described here: [http://software.broadinstitute.org/software/igv/DisplayOptionsAttributes](http://software.broadinstitute.org/software/igv/DisplayOptionsAttributes)

Download the batch script and the attribute file for our dataset:
* Batch script: [Run_batch_IGV_snapshots.txt](http://genomedata.org/gen-viz-workshop/IGV/Run_batch_IGV_snapshots.txt)
* Attribute file: [Igv_HCC1143_attributes.txt](http://genomedata.org/gen-viz-workshop/IGV/Igv_HCC1143_attributes.txt)

Now run the file from the `Tools` menu:

**Automation**
{% include figure.html image="/assets/IGV/run_batch_script.png" width="650" %}

**Notes:**
* This script will navigate automatically to each location in the lab
* A screenshot will be taken and saved to the screenshots directory specified
